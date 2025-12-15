
# Source Code for: A Novel Method for Path-Constrained Dynamic Optimization via Integral Transcription and Inner Approximation

This repository contains the MATLAB implementation for the numerical examples presented in the paper. The code is organized into two main parts corresponding to **Section 5.1** (Single Constraint) and **Section 5.2** (Multiple Constraints).

---

## 📂 Part 1: Single Path Constraint Examples (Section 5.1)

This directory contains the implementation for **Examples 1, 2, and 3**. These examples involve optimal control problems with a single path constraint, comparing the proposed method against several benchmark algorithms.

### 1. Directory Structure

*   **`Proposed_Method/`**: Implementation of the algorithm proposed in this paper (Algorithm 1).
*   **`Comparison_Methods/`**: Implementation of benchmark algorithms for comparison:
    *   `AlphaBB/`, `Interval/`, `Polynomial/`, `Conventional_Discretization/`
*   **`Constraint_Approximation_Comparison/`**: Scripts to reproduce **Figure 7**.

### 2. Prerequisites

*   **MATLAB**: R2024a or later is recommended.
*   **Toolboxes**: Optimization Toolbox (required for `fmincon`).

### 3. How to Run the Examples

#### Mapping of Examples
The code filenames correspond to the paper's examples as follows:

| Paper Section | Code Acronym |
| :--- | :--- |
| **Example 1** | **VDPO** (Van der Pol Oscillator) |
| **Example 2** | **JCB** |
| **Example 3** | **PFBF** (Penicillin Fed-Batch Fermentation) |

#### Execution Modes (Important!)
Each run script (e.g., `run_Ex1_VDPO_Proposed.m`) supports two execution modes, configured at the top of the file:

```matlab
%% --- 1. Configuration ---
ExecutionMode = 'Solve';      % Options: 'Solve', 'Benchmark'
% ExecutionMode = 'Benchmark';
```

*   **`'Solve'` (Default)**: Solves the optimization problem **once**, prints the result, and generates visualization plots.
*   **`'Benchmark'`**: Runs the optimization using MATLAB's `timeit` function to measure accurate CPU time. Use this mode to reproduce the computational performance reported in **Table 1**.
    *   *Note: The `.mat` result files included in the folders contain data from the Benchmark mode.*

#### Running the Code
1.  Navigate to the desired method folder (e.g., `Proposed_Method/`).
2.  Run the corresponding script (e.g., `run_Ex1_VDPO_Proposed.m`).

### 4. Reproducing Figure 7 

To generate **Figure 7**, navigate to the `Constraint_Approximation_Comparison` folder and run:
```matlab
main_plot_constraint_approximations.m
```

### 5. Reproducing Figure 8

Figure 8 in the manuscript illustrates the impact of the number of initial subintervals on the algorithm's performance for Examples 1--3.
To reproduce the **Figure 8**, navigate to the `Proposed_Method` folder and run:

*   `benchmark_grid_performance_Ex1_VDPO.m`
*   `benchmark_grid_performance_Ex2_JCB.m`
*   `benchmark_grid_performance_Ex3_PFBF.m`

---

## 📂 Part 2: Multiple Constraints Example (Section 5.2)

This directory contains the implementation for **Example 4** (Kinematic Car / Parking Problem), which features multiple path constraints. This implementation relies on the **CasADi** framework for automatic differentiation.

### 1. Directory Structure

*   **`run_Parking_Proposed.m`**: The main entry point script to run the example.
*   **`solve_with_UB_CasADi.m`**: The core solver interface linking MATLAB's `fmincon` with CasADi symbolics.

### 2. ⚠️ Critical Prerequisites: CasADi

Unlike the single-constraint examples, this specific example **requires CasADi**. Get the MATLAB binaries from the official website: [https://web.casadi.org/get/](https://web.casadi.org/get/)


### 3. How to Run
1.  Ensure CasADi is added to your path.
2.  Open and run `run_Parking_Proposed.m`.
    *   *Like the previous examples, you can switch between `'Solve'` and `'Benchmark'` modes inside the script.*


---

## 📂 Part 3: Supplementary Tests

This directory contains additional test cases developed to address specific questions raised during the peer review process.

### 1.Oscillatory Constraint Demo (`Example3_Coarse_Initialization_Demo`)

This specific demo was created to verify the algorithm's robustness when handling oscillatory path constraints.

This script forces the proposed algorithm for Example 3 to start with the **coarsest possible partition**: a single time interval covering the entire horizon ($\Pi^0 = \{\Gamma\}$).

#### How to Run
1.  Navigate to the `Example3_Coarse_Initialization_Demo` folder.
2.  Run the script:
    ```matlab
    run_Ex3_Oscillatory_Demo.m
    ```
### 2.Grid Alignment Tests (`Grid_Alignment_Tests_Cases`)

This section contains two test cases (a single-constraint and a multi-constraint problem) designed to demonstrate the algorithm's performance when the discretization of the path constraint is identical to that of the control.

*   **Directories**:
    1.  `03_Supplementary_Tests/Grid_Alignment_Tests_Cases/run_Case1/` (Single Constraint)
    2.  `03_Supplementary_Tests/Grid_Alignment_Tests_Cases/run_Case2/` (Multiple Constraints)

*   **How to Run**:
    Navigate to the respective directory and run the corresponding script (`run_Case1.m` or `run_Case2.m`).

### 3.Discontinuity Handling Demo (`RayleighProblem_Discontinuity_on_Grid`)

A example is included to demonstrate the algorithm's ability to handle path constraints with discontinuous.
Since the control trajectory $u(t)$ is piecewise-constant, the path constraint, $u(t) +\frac{1}{6} x_1(t)$, exhibits discontinuities at the control grid points. 

*   **Directory**: `03_Supplementary_Tests/RayleighProblem_Discontinuity_on_Grid/`
*   **How to Run**: Navigate to the directory and run the script:
    ```matlab
    run_Rayleigh_problem.m
    ```