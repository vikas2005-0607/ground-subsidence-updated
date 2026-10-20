# Updated Ground Subsidence Prediction MATLAB Project

## Files
- `subsidenceGUI.m` — interactive GUI with editable panel start coordinates, calculation extent, and Stop button.
- `subsidence_model.m` — subsidence model updated to honor the user-defined calculation extent and support cooperative cancellation.
- `point_to_line.m` — existing geometry helper.
- `matlab_proj_original.m` — copy of the original script, retained for reference.

## Run
1. Put `subsidenceGUI.m`, `subsidence_model.m`, and `point_to_line.m` in the same MATLAB folder.
2. In MATLAB, change the Current Folder to that folder.
3. Run `subsidenceGUI`.

## New controls
- `Panel start X0/Y0`: lower-left/start coordinate of the panel.
- `Extent beyond panel`: extra distance evaluated outside each of the four panel boundaries, in metres.
- `STOP`: requests cancellation at the next grid-row checkpoint. A stopped run is not drawn or exported as a completed result.

## Contour plots
The 2D subsidence contour uses black contour lines with numeric contour labels and no colour bar. The 3D contour, slope, and strain plots are also line-only rather than filled colour contours.

## Notes
- `Extent beyond panel` changes the calculation grid, not the physical influence-radius equation.
- The Stop button is cooperative: MATLAB processes UI events between grid rows, so stopping occurs at the next checkpoint rather than interrupting a single arithmetic operation.
- These files have been edited programmatically but have not been executed in MATLAB here. Run a test case and compare results with your reference data before using the outputs academically.
