import matplotlib.pyplot as plt
from glob import glob


def get_envelope(file):
    """Get an envelope x, y pairs from an output file
    """
    ts = []
    ps = []
    with open(file) as f:
        for line in f.readlines()[1:]:
            if line.split() == []:
                break
            x, y = line.split()[:2]
            ts.append(float(x))
            ps.append(float(y))
    return ts, ps


# Get all envelout files
files = sorted(glob("bin/envelout*"))

for i, envelout in enumerate(files):
    # Extract relevant data from envelout
    t, p = get_envelope(envelout)
    # Plot data
    plt.plot(t, p, label=f'envelout{i+1}')

plt.legend()
plt.show()
