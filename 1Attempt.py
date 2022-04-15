import os
import numpy as np

Docs = np.empty((5735, 1033))
dic = np.empty((5735,1))
def createMat():
    text_file = open('med.txt', 'r')
    text = text_file.read()
    words = text.split()
    words = [word.strip('.,!;()[]') for word in words]
    words = [word.replace("'s", '') for word in words]
    for elem in words:
        
    print(words)

createMat()