import struct
import numpy as np

def readBinary(fileName):

    with open(fileName, mode='rb') as file: # b is important -> binary
        fileContent = file.read()

    jmax = struct.unpack("i", fileContent[0:4])[0]
    imax = struct.unpack("i", fileContent[4:8])[0]

    print(jmax, imax)

    # data = list( struct.unpack('f'*((len(fileContent) - 24) // 4), fileContent[24:]) )
    # data = np.asarray(data, dtype=np.float32)
    # #data = data.reshape((-1,3))
    # data = data.reshape((-1,2))

    #print(t)
    print(data)

    return data

if __name__ == "__main__":
    fileName = "data.bin"

    data = readBinary(fileName)
