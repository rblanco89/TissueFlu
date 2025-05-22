import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(filename):
    # Load the CSV file
    df = pd.read_csv(filename)

    # Convert time from minutes to days
    df['Days'] = df['Time'] / (60 * 24)

    # Set up the plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Viral Load and IFN
    ax1 = axes[0]
    ax1.set_title('Viral Load and IFN over Time')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Viral Load', color='tab:red')
    ax1.plot(df['Days'], df['ViralLoad'], color='tab:red', label='Viral Load')
    ax1.tick_params(axis='y', labelcolor='tab:red')
    ax1.set_yscale('log')

    ax2 = ax1.twinx()
    ax2.set_ylabel('IFN', color='tab:blue')
    ax2.plot(df['Days'], df['IFN'], color='tab:blue', linestyle='--', label='IFN')
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    ax2.set_yscale('log')

    # Plot 2: Cell States
    axes[1].set_title('Cell States over Time')
    axes[1].set_xlabel('Time (days)')
    axes[1].set_ylabel('Number of Cells')
    axes[1].plot(df['Days'], df['Health'], label='Healthy')
    axes[1].plot(df['Days'], df['Refractory'], label='Refractory')
    axes[1].plot(df['Days'], df['Infected'], label='Infected')
    axes[1].plot(df['Days'], df['Dead'], label='Dead')
    axes[1].legend()
    axes[1].set_yscale('log')

    # Show plot
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python plot_results.py <results_file.csv>")
    else:
        main(sys.argv[1])
