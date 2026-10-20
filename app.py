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
