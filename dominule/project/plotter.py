import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')


def plot_column_histogram(df, column_name):
    plt.figure(figsize=(10, 6))
    plt.hist(df[column_name], bins=30, color='skyblue', edgecolor='black')
    plt.title(f'Histogram of {column_name}')
    plt.xlabel(column_name)
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    plt.show()