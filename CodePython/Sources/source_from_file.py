
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def open_source_file(filepath = "sources/LFXT.txt"):
    headers = ['x', 'y', 'z', 'xs', 'ys', 'E']
    source = pd.read_csv(filepath)
    source.pop(" Charge: -1.5730242e-16 C")
    source = source[1:]
    source.columns = ['DATA']
    
    source[headers] = source.DATA.str.split(" ", expand=True)
    source.pop('DATA')

    data = {}
    for header in headers:
        data[header] = source[header].to_numpy().astype(float)
    return data

def source_from_file(filepath = "sources/LFXT.txt"):
    data = open_source_file(filepath)
    x = data['x']
    y = data['y']
    hist = np.histogram2d(x,y, bins=[20,20])[0]
    remove_columns = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist]
    hist = np.array([i for i, j in zip(hist, remove_columns) if j is False])
    remove_rows = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist.T]
    hist = [i for i, j in zip(hist.T, remove_rows) if j is False]
    return hist

if __name__ == '__main__':
    plt.figure()
    plt.imshow(source_from_file())
    plt.show(block=True)