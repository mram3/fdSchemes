import matplotlib.pyplot as plt
import numpy as np

def extract_data(filename, target_order, target_n):
    """
    Parses the custom CSV format to extract x, exact, numerical, and error lists.
    target_order: "FIRST ORDER" or "SECOND ORDER"
    target_n: 100
    """
    x, exact, numerical, error = [], [], [], []
    capture_flag = False
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Identify the start of the desired data block
        if target_order in line and f"ACCURATE" in line:
            if i + 1 < len(lines) and f"N = {target_n}" in lines[i+1]:
                capture_flag = True
                i += 3  # Skip the 'N=...' line and 'x,exact,numerical,error' header
                continue
                
        # Read data until an empty line is encountered
        if capture_flag:
            if not line:
                break # Reached the end of this data block
            
            vals = line.split(',')
            x.append(float(vals[0]))
            exact.append(float(vals[1]))
            numerical.append(float(vals[2]))
            error.append(float(vals[3]))
            
        i += 1
        
    return np.array(x), np.array(exact), np.array(numerical), np.array(error)

def plot_derivative_analysis(filename, title_prefix):
    # Extract N=100 data for both 1st and 2nd Order schemes
    x_1st, exact_val, num_1st, err_1st = extract_data(filename, "FIRST ORDER", 100)
    x_2nd, _, num_2nd, err_2nd = extract_data(filename, "SECOND ORDER", 100)
    
    if len(x_1st) == 0 or len(x_2nd) == 0:
        print(f"Error: Data for N=100 not found in {filename}.")
        return

    # Create a 1x2 figure for side-by-side comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"{title_prefix} Analysis (N=100)", fontsize=16, fontweight='bold')

    # --- LEFT PLOT: Values ---
    ax1.plot(x_1st, exact_val, 'k-', linewidth=2, label='Exact Analytical')
    ax1.plot(x_1st, num_1st, 'r--', linewidth=1.5, label='1st Order Scheme')
    ax1.plot(x_2nd, num_2nd, 'b:', linewidth=2, label='2nd Order Scheme')
    
    ax1.set_title(f"{title_prefix} - Scheme Values")
    ax1.set_xlabel("x")
    ax1.set_ylabel("Derivative Value")
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='best')

    # --- RIGHT PLOT: Errors (Log Scale) ---
    # We add a small epsilon to the error to avoid log(0) warnings at zero-crossings
    epsilon = 1e-16
    ax2.plot(x_1st, err_1st + epsilon, 'r--', linewidth=1.5, label='1st Order Error')
    ax2.plot(x_2nd, err_2nd + epsilon, 'b:', linewidth=2, label='2nd Order Error')
    
    ax2.set_title(f"{title_prefix} - Absolute Errors (Log Scale)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Absolute Error |Numerical - Exact|")
    ax2.set_yscale('log')
    ax2.grid(True, which="both", linestyle='--', alpha=0.6)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Generate the plots for both output files
    plot_derivative_analysis("First_Derivative.csv", "First Derivative")
    plot_derivative_analysis("Second_Derivative.csv", "Second Derivative")