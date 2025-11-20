# Spline Calculation Tool

This tool is designed to calculate response splines from systematic variations in high-energy physics analysis. It reads particle physics data from ROOT files, applies selections, and for various systematic uncertainties, it calculates splines that model the change in event counts across different kinematic bins.

The output is a ROOT file containing `SplineContainer` objects, which store the binning information and the calculated splines for each systematic.

## Quick Start

1.  **Prepare your data:** You need ROOT files containing the events in a `TTree`.
2.  **Create a configuration file:** Create a `config.yaml` file to specify the input files, binning, systematics, and output file. See the example below.
3.  **Run the tool:**
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ./app/SplineCalcExample ../configs/demo.yaml
    ```

## Configuration

The calculation is controlled by a YAML configuration file. Here is an explanation of the available options with an example:

```yaml
# configs/demo.yaml

# Enable or disable ROOT's implicit multithreading
multithreading: true

# Set the logging level for ROOT's RDataFrame (Debug, Info, Warning, Error)
log_level: "Info"

# Specify the input data files and trees
input_data:
  main_file: "nusyst_genie_extracted.root"
  main_tree: "events"
  friend_files:
    - "nusyst_new_sum.root"
  friend_trees:
    - "sum_tree"

# Define the binning for the analysis
binning_axes:
  - type: "continuous"
    variable: "Ereco_slice"
    edges: [0, 0.5, 1, 1.5, 2, 3, 5, 10]
  - type: "continuous"
    variable: "true_E_from_pdg"
    edges: [0, 1, 2, 3, 4, 5, 10, 20]
  - type: "categorical"
    variable: "isCC"
    categories: [0, 1]

# Set the interpolation type: "Spline" (TSpline3) or "Linear" (TGraph)
interpolation_type: "Spline"

# Information about the systematics
systematics:
  file: "nusyst_genie_extracted.root"
  tree: "meta"

# Name of the output file
output_file: "splines.root"
```

### Configuration Details

-   `multithreading`: If `true`, ROOT's `RDataFrame` will be run in parallel on all available cores.
-   `log_level`: Controls the verbosity of the output from the `RDataFrame`.
-   `input_data`:
    -   `main_file` and `main_tree`: The primary ROOT file and `TTree` containing the event data.
    -   `friend_files` and `friend_trees`: A list of additional files and trees to be "friended" to the main tree. This allows access to their branches as if they were part of the main tree.
-   `binning_axes`: A list of axes that define the histogram bins.
    -   `type`: Can be `continuous` or `categorical`.
    -   `variable`: The name of the branch in the `TTree` to use for this axis.
    -   `edges`: For `continuous` axes, a list of bin edges.
    -   `categories`: For `categorical` axes, a list of integer categories.
-   `interpolation_type`: The type of interpolation to use for the splines. Can be `Spline` for `TSpline3` or `Linear` for `TGraph`.
-   `systematics`:
    -   `file` and `tree`: The file and tree where the systematic parameter information is stored.
-   `output_file`: The path to the output ROOT file where the splines will be saved.

## Code Structure

-   `SplineCalculator`: The main class that orchestrates the entire calculation. It is configured, runs the event loop, and produces the splines.
    -   `SplineCalculator::addBinningAxis()`: Defines a continuous binning axis.
    -   `SplineCalculator::addCategoricalBinningAxis()`: Defines a categorical binning axis.
    -   `SplineCalculator::addSystematic()`: Adds a systematic variation to be calculated.
    -   `SplineCalculator::run()`: Executes the calculation over the input data.
    -   `SplineCalculator::writeSplines()`: Writes the calculated splines to the output file.

-   `SplineContainer`: A container class that holds the results of the calculation for a single systematic. It stores the multidimensional binning (`THnT`) and a map from the global bin number to the calculated spline (`TSpline3` or `TGraph`).

-   `MultiHistoNDAction`: A custom `RDataFrame` action that allows for filling multiple `THnT` (multi-dimensional histograms) in a single, efficient pass over the data. This is used to create one histogram for each systematic parameter variation.

-   `app/SplineCalcExample.cc`: An example executable that shows how to configure and run the `SplineCalculator` using a YAML file.
