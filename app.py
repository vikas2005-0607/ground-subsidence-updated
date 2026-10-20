import io
import math
import textwrap
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(page_title="SubsidenceLab | Ground Subsidence Prediction", page_icon="⛏️", layout="wide")

# ------------------------- Calculation engine -------------------------
RING_WEIGHTS = np.array([0.0297, 0.0854, 0.1311, 0.1607, 0.1708, 0.1596, 0.1294, 0.0859, 0.0397, 0.0069], dtype=float)

def parse_polygon(text):
    """Parse one X,Y vertex per line; accepts comma or whitespace separators."""
    points = []
    for line_no, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.replace(";", ",").split(",")] if "," in line or ";" in line else line.split()
        if len(parts) != 2:
            raise ValueError(f"Polygon line {line_no}: enter exactly X,Y coordinates.")
        try:
            points.append((float(parts[0]), float(parts[1])))
        except ValueError as exc:
            raise ValueError(f"Polygon line {line_no} contains a non-numeric coordinate.") from exc
    if len(points) < 3:
        raise ValueError("A custom polygon needs at least three vertices.")
    arr = np.asarray(points, dtype=float)
    if len(np.unique(arr, axis=0)) < 3:
        raise ValueError("A custom polygon needs at least three distinct vertices.")
    return arr


def point_in_polygon(x, y, polygon):
    """Vectorized ray casting for arrays of points."""
    inside = np.zeros(np.broadcast(x, y).shape, dtype=bool)
    x = np.asarray(x); y = np.asarray(y)
    j = len(polygon) - 1
    for i in range(len(polygon)):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        cross = ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-300) + xi)
        inside ^= cross
        j = i
    return inside


def distance_to_polygon_boundary(x, y, polygon):
    """Minimum Euclidean distance to finite polygon edges, point arrays broadcastable."""
    dmin = np.full(np.broadcast(x, y).shape, np.inf, dtype=float)
    for i in range(len(polygon)):
        ax, ay = polygon[i]
        bx, by = polygon[(i + 1) % len(polygon)]
        vx, vy = bx - ax, by - ay
        denom = vx * vx + vy * vy
        if denom == 0:
            d = np.hypot(x - ax, y - ay)
        else:
            t = np.clip(((x - ax) * vx + (y - ay) * vy) / denom, 0.0, 1.0)
            d = np.hypot(x - (ax + t * vx), y - (ay + t * vy))
        dmin = np.minimum(dmin, d)
    return dmin


def grid_from_bounds(xmin, xmax, ymin, ymax, dx, dy, max_points=220_000):
    xs = np.arange(math.floor(xmin / dx) * dx, math.ceil(xmax / dx) * dx + dx * 0.25, dx)
    ys = np.arange(math.floor(ymin / dy) * dy, math.ceil(ymax / dy) * dy + dy * 0.25, dy)
    if len(xs) * len(ys) > max_points:
        raise ValueError(f"Grid has {len(xs)*len(ys):,} points. Increase dx/dy or reduce extent (maximum {max_points:,} points).")
    return np.meshgrid(xs, ys)


def generalized_model(p):
    """Vectorized port of subsidence_model.m, preserving its 10-ring/64-sector equations."""
    H, W, L = p["H"], p["W"], p["L"]
    draw, NEW, h, e = p["drawDeg"], p["NEW"], p["h"], p["e"]
    x0, y0, dx, dy = p["x0"], p["y0"], p["dx"], p["dy"]
    mode, T = p["mode"], p["T"]
    panel = p["polygon"]
    xminp, yminp = panel.min(axis=0); xmaxp, ymaxp = panel.max(axis=0)
    panel_width, panel_length = xmaxp - xminp, ymaxp - yminp
    if H <= 0 or W <= 0 or L <= 0 or NEW <= 0 or h <= 0 or e <= 0 or draw <= 0:
        raise ValueError("H, W, L, NEW, h, e and draw angle must be positive.")
    if e > 1: raise ValueError("Extraction factor e must be between 0 and 1.")
    if draw >= 90: raise ValueError("Draw angle must be less than 90 degrees.")
    if mode == "Multi-seam" and not (0 <= T < H): raise ValueError("For multi-seam mode, 0 ≤ T < H is required.")
    if panel_width <= 0 or panel_length <= 0: raise ValueError("Panel polygon must have non-zero width and length.")

    R = H * math.tan(math.radians(draw))
    ext = max(p["extent"], 0.0)
    # The supplied GUI passes extent, but its MATLAB model did not consume it. This port makes it effective.
    pad = max(R, ext)
    xmin = max(0.0, xminp - pad); xmax = xmaxp + pad
    ymin = max(0.0, yminp - pad); ymax = ymaxp + pad
    X, Y = grid_from_bounds(xmin, xmax, ymin, ymax, dx, dy)
    flatx, flaty = X.ravel(), Y.ravel()
    n = flatx.size
    S = np.zeros(n); Ux = np.zeros(n); Uy = np.zeros(n)

    width_ratio = panel_width / (H * NEW)
    if mode == "Single Seam":
        base_factor = 0.33 * (1 + 1.1 * math.tanh(1.4 * (width_ratio - 1.8)))
    else:
        base_factor = 0.12 * (1 + 5.9 * math.tanh(0.38 * (width_ratio - 1.45)))
    base_factor = max(base_factor, 0.0)
    Smax = base_factor * h * e
    sector_weights = RING_WEIGHTS / 64.0
    dr = R / 10.0
    rho = (np.arange(1, 11) - 0.5) * dr
    theta = (np.arange(1, 65) - 0.5) * (2 * np.pi / 64.0)
    ct, stheta = np.cos(theta), np.sin(theta)

    # Chunking bounds memory for high-resolution grids.
    chunk = 1800
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        px = flatx[start:stop, None]; py = flaty[start:stop, None]
        sz = np.zeros(stop - start); sx = np.zeros(stop - start); sy = np.zeros(stop - start)
        for ri, radius in enumerate(rho):
            ex = px + radius * ct[None, :]
            ey = py + radius * stheta[None, :]
            inside = point_in_polygon(ex, ey, panel)
            if not inside.any():
                continue
            db = distance_to_polygon_boundary(ex, ey, panel)
            # Preserve the MATLAB GUI's distance from its starting face and dmax convention.
            d = ey - yminp
            dmax = panel_length
            local_new = np.full(ex.shape, NEW)
            local_smax = np.full(ex.shape, Smax)
            if mode == "Multi-seam":
                gx0, gy0 = p["goafX0"], p["goafY0"]
                gx1, gy1 = gx0 + p["goafW"], gy0 + p["goafL"]
                ingoaf = (ex >= gx0) & (ex <= gx1) & (ey >= gy0) & (ey <= gy1)
                local_new[ingoaf] = (0.6 * NEW * (H - T) + NEW * T) / H
                local_smax[ingoaf] = 1.17 * base_factor * h * e
            wz = 0.5 * np.tanh((5 * db) / (1.5 * local_new * h) - 2.4) + 0.5
            qz = 0.9 - 0.1 * np.tanh(0.5 * (d - 0.4 * dmax))
            dvz = -sector_weights[ri] * local_smax * wz * qz
            dvz = np.where(inside, dvz, 0.0)
            sz += dvz.sum(axis=1)
            sx += (dvz * (ex - px) / H).sum(axis=1)
            sy += (dvz * (ey - py) / H).sum(axis=1)
        S[start:stop] = sz; Ux[start:stop] = sx; Uy[start:stop] = sy

    S = S.reshape(X.shape); Ux = Ux.reshape(X.shape); Uy = Uy.reshape(X.shape)
    # MATLAB gradient(F,dy,dx) returns derivatives along row/y and column/x.
    dSdy, dSdx = np.gradient(S, dy, dx)
    slope = np.hypot(dSdx, dSdy)
    dUxdy, dUxdx = np.gradient(Ux, dy, dx)
    dUydy, dUydx = np.gradient(Uy, dy, dx)
    Exx = dUxdx; Eyy = dUydy; Exy = 0.5 * (dUxdy + dUydx)
    mean_strain = 0.5 * (Exx + Eyy)
    radius_strain = np.sqrt(((Exx - Eyy) / 2) ** 2 + Exy ** 2)
    principal1 = mean_strain + radius_strain
    principal2 = mean_strain - radius_strain
    if width_ratio < 1: width_class = "Non-effective"
    elif width_ratio < 3.5: width_class = "Sub-critical"
    else: width_class = "Super-critical"
    return dict(X=X, Y=Y, Subsidence=S, Ux=Ux, Uy=Uy, slopeMagnitude=slope,
                principal1=principal1, principal2=principal2, Exx=Exx, Eyy=Eyy, Exy=Exy,
                R=R, widthRatio=width_ratio, widthClass=width_class, baseFactor=base_factor,
                Smax=Smax, panel=panel, maxTensileStrain=float(np.max(principal1)),
                maxCompressiveStrain=float(np.min(principal2)), dx=dx, dy=dy,
                method="GUI model (ported from subsidence_model.m)")


def original_script_model(p):
    """Port of matlab_proj.m legacy script: fixed default panel/grid/sampling and its original scale factor."""
    H, W, draw, NEW, h, e = p["H"], p["W"], p["drawDeg"], p["NEW"], p["h"], p["e"]
    if min(H, W, draw, NEW, h, e) <= 0: raise ValueError("H, W, draw angle, NEW, h and e must be positive.")
    if draw >= 90: raise ValueError("Draw angle must be less than 90 degrees.")
    if e > 1: raise ValueError("Extraction factor e must be between 0 and 1.")
    # Original file hard-codes the panel coordinates, grid, and dmax; preserve this path intentionally.
    x0, y0, panelW, panelL = 100.0, 100.0, 198.0, 294.0
    panel = np.array([[100,100],[298,100],[298,394],[100,394]], dtype=float)
    xs = np.arange(1.0, 401.0, 4.0); ys = np.arange(1.0, 497.0, 5.0)
    X, Y = np.meshgrid(xs, ys)
    flatx, flaty = X.ravel(), Y.ravel(); n = flatx.size
    S = np.zeros(n)
    R = H * math.tan(math.radians(draw)); r = R / 10.0
    sub_factor = 0.33 * (1 + 1.11 * math.tanh(1.4 * (W / H / NEW - 1.8)))
    Smax = sub_factor * h * e
    ring_weights = RING_WEIGHTS / 65.0
    theta = np.arange(0, 65) * (np.pi / 32.0)
    ct, stheta = np.cos(theta), np.sin(theta)
    chunk = 1600
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        px = flatx[start:stop, None]; py = flaty[start:stop, None]
        total = np.zeros(stop - start)
        for i in range(1, 11):
            radius = i * r
            ex = px + radius * ct[None, :]
            ey = py + radius * stheta[None, :]
            inside = point_in_polygon(ex, ey, panel)
            if not inside.any(): continue
            # Original uses perpendicular distances to the rectangle's four infinite edge lines.
            db = np.minimum.reduce([np.abs(ex - 100.0), np.abs(ex - 298.0), np.abs(ey - 100.0), np.abs(ey - 394.0)])
            d = ey - 100.0
            dmax = 354.46
            wz = 0.5 * np.tanh((5 * db) / (1.5 * NEW * h) - 2.4) + 0.5
            qz = 0.9 - 0.1 * np.tanh(0.5 * (d - 0.4 * dmax))
            contrib = ring_weights[i-1] * (-Smax) * wz * qz * 100.0
            total += np.where(inside, contrib, 0.0).sum(axis=1)
        S[start:stop] = total
    S = S.reshape(X.shape)
    dSdy, dSdx = np.gradient(S, 5.0, 4.0)
    slope = np.hypot(dSdx, dSdy)
    return dict(X=X, Y=Y, Subsidence=S, slopeMagnitude=slope, R=R,
                widthRatio=W/(H*NEW), widthClass="Legacy standalone script",
                Smax=Smax, panel=panel, dx=4.0, dy=5.0,
                method="Original standalone script (ported from matlab_proj.m)")

# ------------------------------ UI -------------------------------------
def apply_style():
    st.markdown("""
    <style>
    .stApp {background:#f6f8fb;color:#182230;}
    [data-testid="stHeader"] {background:rgba(246,248,251,.92);}
    .block-container {padding-top:1.4rem;padding-bottom:2.5rem;max-width:1500px;}
    .hero {background:linear-gradient(115deg,#0d2748,#164e83);padding:1.45rem 1.7rem;border-radius:18px;color:white;margin-bottom:1rem;}
    .hero h1 {font-size:2rem;margin:0 0 .3rem 0;color:white;}
    .hero p {margin:0;color:#d8e8f8;}
    .section-card {background:white;border:1px solid #e1e7ef;border-radius:14px;padding:1rem 1.1rem;margin:.25rem 0 1rem 0;}
    .step-title {font-size:1.08rem;font-weight:750;color:#12365b;margin:.15rem 0 .7rem 0;}
    div.stButton > button[kind="primary"] {background:#1264a3;border-color:#1264a3;color:white;font-weight:700;}
    div[data-testid="stMetric"] {background:white;border:1px solid #e4eaf2;padding:.8rem;border-radius:12px;}
    </style>""", unsafe_allow_html=True)

apply_style()

with st.sidebar:
    st.markdown("## ⛏️ SubsidenceLab")
    page = st.radio("Navigation", ["Predict", "About & Method"], label_visibility="collapsed")
    st.markdown("---")
    st.caption("Ground Subsidence Prediction • Mining Engineering")
    st.caption("Academic model: outputs require verification against field measurements.")

st.markdown('<div class="hero"><h1>SubsidenceLab</h1><p>Ground Subsidence Prediction · Modified influence-function method</p></div>', unsafe_allow_html=True)

if page == "About & Method":
    st.markdown("### About the model")
    st.write("This application ports the supplied MATLAB workflows into Python for browser-based use. The GUI method and legacy standalone script are kept as separate calculation paths because the source implementations use different ring sampling and scaling.")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("#### GUI model — `subsidence_model.m`")
        st.markdown("- 10 concentric rings and 64 sector-centroid samples per ring.\n- Uses the supplied ring weights normalized by 64.\n- Computes subsidence, horizontal displacement, slope and numerical principal strains.\n- Supports single-seam and multi-seam correction as encoded in the MATLAB source.")
    with c2:
        st.markdown("#### Original script — `matlab_proj.m`")
        st.markdown("- Uses the original fixed 100 × 100 grid.\n- Uses 10 radii and 65 angular samples per ring.\n- Preserves the original `*100` contribution scale and fixed panel geometry.\n- Its outputs are not expected to match the GUI model exactly.")
    st.warning("Validation status: the source files were translated, but the Python results have not been executed and compared with MATLAB outputs in this environment. Verify identical test cases before academic or engineering use.")
    st.markdown("#### Model assumptions and limitations")
    st.write("The GUI source uses a rectangular extraction panel. This app includes a custom-polygon option as an extension using point-in-polygon tests and distance to finite polygon edges; that extension is not part of the original MATLAB implementation and should be independently validated. The GUI model's extent setting controls the output grid padding; the source MATLAB model passed this parameter but did not consume it in the grid bounds, so this app makes it functional.")
else:
    # input controls
    st.markdown('<div class="section-card"><div class="step-title">01 · Calculation method</div>', unsafe_allow_html=True)
    method = st.radio("Select calculation method", ["GUI model", "Original MATLAB script"], horizontal=True,
                      help="The two MATLAB source paths use different sampling and scaling. Results may differ.")
    st.markdown('</div>', unsafe_allow_html=True)

    left, right = st.columns([1.02, 1.4], gap="large")
    with left:
        st.markdown('<div class="section-card"><div class="step-title">02 · Mining panel geometry</div>', unsafe_allow_html=True)
        geom_mode = st.radio("Panel geometry", ["Rectangle", "Custom polygon"], horizontal=True,
                             disabled=(method == "Original MATLAB script"),
                             help="The original script has a hard-coded rectangular panel. Custom polygon is an extension for the GUI model.")
        if geom_mode == "Rectangle" or method == "Original MATLAB script":
            c1, c2 = st.columns(2)
            with c1:
                x0 = st.number_input("Panel start X₀ (m)", value=100.0, step=10.0)
                W = st.number_input("Panel width W (m)", min_value=1.0, value=198.0, step=10.0)
            with c2:
                y0 = st.number_input("Panel start Y₀ (m)", value=100.0, step=10.0)
                L = st.number_input("Panel length L (m)", min_value=1.0, value=294.0, step=10.0)
            polygon = np.array([[x0,y0],[x0+W,y0],[x0+W,y0+L],[x0,y0+L]], dtype=float)
            polygon_text = ""
        else:
            st.caption("Enter one vertex per line as X,Y. Vertices should be ordered around the panel boundary.")
            polygon_text = st.text_area("Panel vertices (m)", value="100,100\n298,100\n298,394\n100,394", height=110)
            try:
                polygon = parse_polygon(polygon_text)
                xmin_poly, ymin_poly = polygon.min(axis=0); xmax_poly, ymax_poly = polygon.max(axis=0)
                x0, y0, W, L = float(xmin_poly), float(ymin_poly), float(xmax_poly-xmin_poly), float(ymax_poly-ymin_poly)
            except ValueError as err:
                st.error(str(err)); polygon = np.array([[0,0],[1,0],[1,1],[0,1]], float)
                x0=y0=0.0; W=L=1.0
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-card"><div class="step-title">03 · Model parameters</div>', unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        with c1:
            H = st.number_input("Depth of cover H (m)", min_value=0.1, value=300.0, step=10.0)
            drawDeg = st.number_input("Angle of draw γ (°)", min_value=0.1, max_value=89.0, value=13.0, step=0.5)
            NEW = st.number_input("NEW parameter", min_value=0.001, value=0.40, step=0.05, format="%.3f")
        with c2:
            h = st.number_input("Seam thickness h (m)", min_value=0.1, value=4.6, step=0.1)
            e = st.number_input("Extraction factor e", min_value=0.001, max_value=1.0, value=0.66, step=0.01)
            mode = st.selectbox("Mining mode", ["Single Seam", "Multi-seam"], disabled=(method == "Original MATLAB script"))
        c3, c4 = st.columns(2)
        with c3:
            T = st.number_input("Parting T (m)", min_value=0.0, value=50.0, step=5.0, disabled=(mode != "Multi-seam" or method == "Original MATLAB script"))
            goafW = st.number_input("Goaf width (m)", min_value=1.0, value=float(W), step=10.0, disabled=(mode != "Multi-seam" or method == "Original MATLAB script"))
        with c4:
            goafL = st.number_input("Goaf length (m)", min_value=1.0, value=float(L), step=10.0, disabled=(mode != "Multi-seam" or method == "Original MATLAB script"))
            st.caption("Influence radius R = H × tan(γ)")
            st.metric("Influence radius R", f"{H * math.tan(math.radians(drawDeg)):.2f} m")
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-card"><div class="step-title">04 · Calculation boundary & grid</div>', unsafe_allow_html=True)
        if method == "Original MATLAB script":
            st.info("Legacy script retains its fixed X=1:4:400 and Y=1:5:496 grid and fixed panel coordinates, as in the original source.")
            extent = 100.0; dx = 4.0; dy = 5.0
        else:
            c1, c2, c3 = st.columns(3)
            with c1: extent = st.number_input("Extent beyond panel (m)", min_value=0.0, value=100.0, step=10.0)
            with c2: dx = st.number_input("Grid spacing dx (m)", min_value=1.0, value=4.0, step=1.0)
            with c3: dy = st.number_input("Grid spacing dy (m)", min_value=1.0, value=5.0, step=1.0)
            st.caption("Boundary coordinates can be pasted into the polygon editor above. Rectangle mode generates the boundary automatically.")
        st.markdown('</div>', unsafe_allow_html=True)
        b1, b2 = st.columns(2)
        with b1: calculate = st.button("⛏️  CALCULATE", type="primary", use_container_width=True)
        with b2: reset = st.button("↺  RESET", use_container_width=True)
        if reset:
            st.session_state.clear()
            st.rerun()

    with right:
        st.markdown('<div class="section-card"><div class="step-title">05 · Prediction results</div>', unsafe_allow_html=True)
        if calculate:
            try:
                if method == "Original MATLAB script":
                    params = dict(H=H, W=198.0, drawDeg=drawDeg, NEW=NEW, h=h, e=e)
                    result = original_script_model(params)
                else:
                    if geom_mode == "Custom polygon":
                        polygon = parse_polygon(polygon_text)
                    params = dict(H=H, W=W, L=L, drawDeg=drawDeg, NEW=NEW, h=h, e=e,
                                  x0=x0, y0=y0, dx=dx, dy=dy, mode=mode, T=T,
                                  extent=extent, polygon=polygon,
                                  goafX0=x0, goafY0=y0, goafW=goafW, goafL=goafL)
                    result = generalized_model(params)
                st.session_state["subsidence_result"] = result
                st.session_state["subsidence_method"] = method
                st.success("Calculation completed. Review the results and validate against MATLAB before research use.")
            except Exception as exc:
                st.error(f"Calculation failed: {exc}")
        result = st.session_state.get("subsidence_result")
        if result is None:
            st.info("Set the inputs on the left and click CALCULATE to generate model results.")
        else:
            z = result["Subsidence"]
            k1, k2, k3, k4 = st.columns(4)
            k1.metric("Maximum downward magnitude", f"{np.max(np.abs(z)):.5f} m")
            k2.metric("Minimum subsidence value", f"{np.min(z):.5f} m")
            k3.metric("Influence radius", f"{result['R']:.2f} m")
            k4.metric("Grid points", f"{z.size:,}")
            st.caption(f"Method: {result['method']} · Width classification: {result['widthClass']} · Width ratio: {result['widthRatio']:.3f}")
        st.markdown('</div>', unsafe_allow_html=True)

    result = st.session_state.get("subsidence_result")
    if result is not None:
        st.markdown("### Results visualization")
        tabs = st.tabs(["2D Contour", "3D Surface", "Profiles", "Displacement / Slope", "Strain", "Data & Export"])
        X, Y, Z = result["X"], result["Y"], result["Subsidence"]
        with tabs[0]:
            fig = go.Figure()
            fig.add_trace(go.Contour(x=X[0,:], y=Y[:,0], z=Z, contours=dict(coloring="lines", showlabels=True),
                                     line=dict(color="#111827", width=1), showscale=False, name="Subsidence"))
            panel = result["panel"]
            fig.add_trace(go.Scatter(x=np.r_[panel[:,0],panel[0,0]], y=np.r_[panel[:,1],panel[0,1]], mode="lines",
                                     line=dict(color="#dc2626", width=2, dash="dash"), name="Mining panel"))
            fig.update_layout(template="plotly_white", height=590, xaxis_title="X (m)", yaxis_title="Y (m)",
                              title="Ground subsidence — labelled contour lines", margin=dict(l=20,r=20,t=55,b=20))
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            st.plotly_chart(fig, use_container_width=True)
        with tabs[1]:
            fig3 = go.Figure(data=[go.Surface(x=X, y=Y, z=Z, colorscale="Blues", colorbar=dict(title="m"), contours={
                "z": {"show": True, "usecolormap": False, "highlightcolor": "#111827", "project": {"z": True}}})])
            panel = result["panel"]
            fig3.add_trace(go.Scatter3d(x=np.r_[panel[:,0],panel[0,0]], y=np.r_[panel[:,1],panel[0,1]], z=np.zeros(len(panel)+1),
                                       mode="lines", line=dict(color="red", width=5), name="Mining panel"))
            fig3.update_layout(template="plotly_white", height=650, scene=dict(xaxis_title="X (m)", yaxis_title="Y (m)", zaxis_title="Subsidence (m)"),
                               title="3D subsidence surface")
            st.plotly_chart(fig3, use_container_width=True)
        with tabs[2]:
            c1, c2 = st.columns(2)
            with c1:
                row_index = st.slider("Select Y row", 0, len(Y[:,0])-1, len(Y[:,0])//2, key="profile_y")
                pf = go.Figure(go.Scatter(x=X[row_index,:], y=Z[row_index,:], mode="lines", line=dict(color="#1264a3", width=2)))
                pf.update_layout(template="plotly_white", title=f"X profile at Y = {Y[row_index,0]:.2f} m", xaxis_title="X (m)", yaxis_title="Subsidence (m)", height=380)
                st.plotly_chart(pf, use_container_width=True)
            with c2:
                col_index = st.slider("Select X column", 0, len(X[0,:])-1, len(X[0,:])//2, key="profile_x")
                pf2 = go.Figure(go.Scatter(x=Y[:,col_index], y=Z[:,col_index], mode="lines", line=dict(color="#0f766e", width=2)))
                pf2.update_layout(template="plotly_white", title=f"Y profile at X = {X[0,col_index]:.2f} m", xaxis_title="Y (m)", yaxis_title="Subsidence (m)", height=380)
                st.plotly_chart(pf2, use_container_width=True)
        with tabs[3]:
            if "Ux" in result:
                c1, c2 = st.columns(2)
                with c1:
                    figu = go.Figure(go.Contour(x=X[0,:], y=Y[:,0], z=result["Ux"], contours=dict(coloring="lines", showlabels=True), line=dict(color="#111827"), showscale=False))
                    figu.update_layout(template="plotly_white", title="Horizontal displacement Ux (m)", height=450)
                    st.plotly_chart(figu, use_container_width=True)
                with c2:
                    figs = go.Figure(go.Contour(x=X[0,:], y=Y[:,0], z=result["slopeMagnitude"], contours=dict(coloring="lines", showlabels=True), line=dict(color="#111827"), showscale=False))
                    figs.update_layout(template="plotly_white", title="Slope magnitude", height=450)
                    st.plotly_chart(figs, use_container_width=True)
            else:
                st.info("The original standalone script does not compute horizontal displacement in the supplied source.")
        with tabs[4]:
            if "principal1" in result:
                c1, c2 = st.columns(2)
                with c1:
                    f1 = go.Figure(go.Contour(x=X[0,:], y=Y[:,0], z=result["principal1"], contours=dict(coloring="lines", showlabels=True), line=dict(color="#111827"), showscale=False))
                    f1.update_layout(template="plotly_white", title="Maximum principal strain", height=450)
                    st.plotly_chart(f1, use_container_width=True)
                with c2:
                    f2 = go.Figure(go.Contour(x=X[0,:], y=Y[:,0], z=result["principal2"], contours=dict(coloring="lines", showlabels=True), line=dict(color="#111827"), showscale=False))
                    f2.update_layout(template="plotly_white", title="Minimum principal strain", height=450)
                    st.plotly_chart(f2, use_container_width=True)
                st.write(f"Maximum principal strain: **{result['maxTensileStrain']:.6g}** · Minimum principal strain: **{result['maxCompressiveStrain']:.6g}**")
            else:
                st.info("Principal strain is not calculated by the original standalone script.")
        with tabs[5]:
            df = pd.DataFrame({"X_m": X.ravel(), "Y_m": Y.ravel(), "Subsidence_m": Z.ravel()})
            for key, label in [("Ux","Ux_m"),("Uy","Uy_m"),("slopeMagnitude","Slope"),("principal1","Principal1"),("principal2","Principal2")]:
                if key in result: df[label] = result[key].ravel()
            st.dataframe(df.head(1000), use_container_width=True, height=350)
            st.caption(f"Showing first {min(1000, len(df)):,} rows of {len(df):,}. Download includes all grid points.")
            st.download_button("Download all results (CSV)", df.to_csv(index=False).encode("utf-8"), file_name="subsidence_results.csv", mime="text/csv")
            if "Ux" in result:
                st.write(f"Max tensile strain: {result['maxTensileStrain']:.6g}; max compressive strain: {result['maxCompressiveStrain']:.6g}")

st.markdown("---")
st.caption("SubsidenceLab · Research prototype. Verify units, model equations, and predictions against original MATLAB output and observed data before engineering decisions.")
