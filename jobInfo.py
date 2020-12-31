def getJobParams(outputDir):

    if outputDir == 3:
        return [1e-4, 10, 30]

    # Dec 14
    alphaL = [1e-3, 1e-4]
    amplitudeL = [10, 30]
    positionL = [30, 60, 90]
    startingDir = 30
    for i in alphaL:
        for j in amplitudeL:
            for k in positionL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    position = k
                    return [alpha, amplitude, position]
                else:
                    startingDir += 1

    # Dec 17  - non inverted
    startingDir = 45
    amplitudeL = [30]
    positionL = [10, 20, 30, 60, 90]
    for i in alphaL:
        for j in amplitudeL:
            for k in positionL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    position = k
                    return [alpha, amplitude, position]
                else:
                    startingDir += 1

    # Dec 17 - inverted
    startingDir = 55
    amplitudeL = [3, 30]
    positionL = [50, 100]
    for i in alphaL:
        for j in amplitudeL:
            for k in positionL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    position = k
                    return [alpha, amplitude, position]
                else:
                    startingDir += 1

    # Dec 23 stationary scripts
    alphaL = [1e-3, 1e-4]
    amplitudeL = [10, 30]
    positionL = [30, 60, 90]
    startingDir = 70
    for i in alphaL:
        for j in amplitudeL:
            for k in positionL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    position = k
                    return [alpha, amplitude, position]
                else:
                    startingDir += 1

    return [0,0,0]
