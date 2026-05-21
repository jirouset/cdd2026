import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')


def plot_column_histogram(df, column_name, threshold=None):
    plt.figure(figsize=(10, 6))
    plt.hist(df[column_name], bins=30, color='skyblue', edgecolor='black')

    if threshold:
        plt.axvline(threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold: {threshold} nM')

    plt.title(f'Histogram of {column_name}')
    plt.xlabel(column_name)
    plt.ylabel('Frequency')
    plt.legend()  # Displays the threshold label
    plt.grid(axis='y', alpha=0.75)
    plt.show()