import matplotlib.pyplot as plt
import csv

#Initialising lists to store data from csv
x = []
numerical = []
exact = []
error = []
plot_title =""

#Reading the csv file
with open('1_PolynomialDS_Output.csv', 'r') as file:
    reader = csv.reader(file)

    title_row = next(reader)
    if title_row:
        plot_title = title_row[0]

    next(reader) #skipping the column heading row

    #reading the data in csv file
    for row in reader:
        if not row: continue #skipping any empty line
        x.append(float(row[0]))
        numerical.append(float(row[1]))
        exact.append(float(row[2]))
        error.append(float(row[3]))

#setting up the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,5))
fig.suptitle(plot_title, fontsize=16, fontweight='bold')

ax1.plot(x, exact, label='Exact Solution', color ='tab:orange', linewidth=2.5)
ax1.plot(x, numerical, label='Numerical Solution', color='tab:blue', marker='o', linestyle='--', markersize=4)

ax1.set_title('Numerical vs Exact Solution', fontsize=12)
ax1.set_xlabel('x', fontsize = 11)
ax1.set_ylabel('Derivative Value', fontsize=11)
ax1.grid(True, linestyle=':',alpha=0.7)
ax1.legend()

ax2.plot(x, error, label='Absolute Error', color='tab:blue', marker='o', linestyle='-', markersize=4)

ax2.set_title('Truncation Error Profile', fontsize=12)
ax2.set_xlabel('x', fontsize=11)
ax2.set_ylabel('Absolute Error', fontsize=11)
ax2.grid(True, linestyle=':', alpha=0.7)
ax2.legend()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.tight_layout()
plt.savefig("1_PolynomialDS_plot.png", dpi = 300)
plt.show()