# M5 Experimental Design Assistant for Bxb1 Integration

A predictive model for optimizing Bxb1 serine integrase-mediated genomic integration.

## Quick Start

### Installation

git clone https://github.com/[NYCU-Formosa]/M5-Bxb1-Design-Tool.git
cd M5-Bxb1-Design-Tool
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt

### Basic Usage

from m5_tool import M5Model

# Load model
model = M5Model()

# Predict success rate for 1.5 kb insert
P = model.predict_mode_b(L=1.5)  # L in kb
print(f"Success rate: {P*100:.1f}%")

# Calculate expected colonies
colonies = model.expected_colonies(L=1.5, N_CFU=30)
print(f"Expected colonies: {colonies:.1f} per plate")

### Generate Figures

# Generate specific figure
python plot_figures.py --figure 1 --output fig1.png

# Generate all figures
python plot_figures.py --figure 0


### Batch Query

python m5_tool.py --input data/example_queries.csv --output results.csv


## Calibration for Your Strain

from m5_tool import M5Model

model = M5Model()
mae = model.calibrate('your_data.csv', params_to_fit=['a', 'b'])
print(f"MAE: {mae:.2f}%")
model.save('custom_params.json')


## Data Format

### Input CSV for calibration

Your CSV file should have these columns:

- `insert_length_bp`: Insert size in base pairs
- `success_rate`: Observed integration success rate (0-1)
- `n_replicates`: Number of biological replicates

Example:

insert_length_bp,success_rate,n_replicates
1039,0.875,3
1605,0.790,3
1613,0.800,3


## Model Parameters

The model uses the following calibrated parameters:

- **a** = 3.0161 (logit intercept)
- **b** = 1.03×10⁻³ /bp (length sensitivity)
- **K** = 0.20 μM (concentration half-saturation)
- **C_ref** = 0.387 μM (reference concentration)

These parameters were calibrated on our strain. You can recalibrate for your own strain using the calibration script.

## Citation

If you use this tool in your research, please cite:

[NYCU-Formosa]. (2025). M5 Experimental Design Assistant for Bxb1 Integration. 
iGEM 2025. https://github.com/[YourTeam]/M5-Bxb1-Design-Tool

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open an issue on GitHub or contact [nycuformosa.igem@gmail.com]

## Acknowledgments

Developed as part of iGEM 2025 project.
