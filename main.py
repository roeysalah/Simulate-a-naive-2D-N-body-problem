import matplotlib.pyplot as plt

def scatter_data_from_files(file,title,save_file):
    # Read data from the text file
    data = []
    with open(file, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split())
            data.append((x, y))

    # Scatter plot the data
    plt.scatter(*zip(*data))

    # Add labels and legend
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

    # Set the title
    plt.title(title)
    plt.savefig(save_file)

    # Show the plot
    plt.show()

# Scatter plot for start positions
scatter_data_from_files('snapshot_t0.txt', 'Galaxy at T = 0 ','Galaxy0.png')

# Scatter plot for middle positions
scatter_data_from_files('snapshot_t5000.txt', 'Galaxy at T = Tmax / 2 ','GalaxyhalfTmax.png')

# Scatter plot for end positions
scatter_data_from_files('snapshot_t10000.txt',  'Galaxy at T = Tmax ','GalaxyTmax.png')
