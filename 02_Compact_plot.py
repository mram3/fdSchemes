import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_results(file_path):
    # Initialize lists to hold the parsed data
    first_deriv_data = []
    second_deriv_data = []
    current_section = None
    headers = []

    # Read and parse the custom CSV format
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            # Detect section changes
            if line.startswith("First Derivative"):
                current_section = 1
                continue
            elif line.startswith("Second Derivative"):
                current_section = 2
                continue
            elif line.startswith("Performance"):
                continue
                
            # Capture headers or data
            if line.startswith("x,"):
                headers = line.split(',')
            else:
                # Convert comma-separated string values to floats
                values = [float(v) for v in line.split(',')]
                if current_section == 1:
                    first_deriv_data.append(values)
                elif current_section == 2:
                    second_deriv_data.append(values)

    # Convert to Pandas DataFrames for easy manipulation
    df1 = pd.DataFrame(first_deriv_data, columns=headers)
    df2 = pd.DataFrame(second_deriv_data, columns=headers)

    return df1, df2

def calculate_errors(df):
    # Calculate absolute errors compared to Exact solution
    df['Pade Error'] = np.abs(df['Pade4'] - df['Exact'])
    df['FD Error'] = np.abs(df['FourthOrder FD'] - df['Exact'])
    
    # Avoid zero values for log plots by substituting with machine epsilon
    epsilon = np.finfo(float).eps
    df['Pade Error'] = df['Pade Error'].replace(0, epsilon)
    df['FD Error'] = df['FD Error'].replace(0, epsilon)
    
    return df

def plot_side_by_side(df, title):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Left Subplot: Values ----
    ax1.plot(df['x'], df['Exact'], 'k-', label='Exact', linewidth=2)
    ax1.plot(df['x'], df['Pade4'], 'r--', label='Compact (Pade4)', linewidth=1.5)
    ax1.plot(df['x'], df['FourthOrder FD'], 'b:', label='Conventional FD 4th', linewidth=1.5)
    
    ax1.set_title(f"{title} - Scheme Values")
    ax1.set_xlabel("x")
    ax1.set_ylabel("Derivative Value")
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.6)

    # ---- Right Subplot: Errors ----
    ax2.plot(df['x'], df['Pade Error'], 'r--', label='Pade4 Error', linewidth=1.5)
    ax2.plot(df['x'], df['FD Error'], 'b:', label='FD Error', linewidth=1.5)
    
    # Log scale is usually best for viewing numerical scheme errors
    ax2.set_yscale('log') 
    
    ax2.set_title(f"{title} - Absolute Errors (Log Scale)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Absolute Error |Numerical - Exact|")
    ax2.legend()
    ax2.grid(True, which="both", linestyle='--', alpha=0.6)

    plt.suptitle(title, fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # 1. Parse the generated CSV file
    file_name = "02_Compact_results.csv"
    try:
        df_first, df_second = parse_results(file_name)
    except FileNotFoundError:
        print(f"Error: Could not find '{file_name}'. Ensure you have run the C++ program first.")
        exit()

    # 2. Calculate errors for both derivatives
    df_first = calculate_errors(df_first)
    df_second = calculate_errors(df_second)

    # Print a quick tabular summary
    print("\n--- First Derivative Sample Data ---")
    print(df_first[['x', 'Exact', 'Pade4', 'FourthOrder FD', 'Pade Error', 'FD Error']].head())

    print("\n--- Second Derivative Sample Data ---")
    print(df_second[['x', 'Exact', 'Pade4', 'FourthOrder FD', 'Pade Error', 'FD Error']].head())

    # 3. Generate side-by-side plots for First and Second Derivatives
    plot_side_by_side(df_first, "First Derivative Analysis")
    plot_side_by_side(df_second, "Second Derivative Analysis")