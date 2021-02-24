def getJobParams(outputDir):

    # wide
    if outputDir == 232:
        return [1e-3, 10, 90]
    elif outputDir == 233 or outputDir == 234:
        return [1e-3, 10, 120]

    # vfrag
    alphaL = [1e-4]
    amplitudeL = [10, 30]
    posL = [30, 60, 90]
    startingDir = 220
    for i in alphaL:
        for j in amplitudeL:
            for k in posL:
                for vfrag in amplitudeL:
                    if startingDir == outputDir:
                        alpha = i
                        amplitude = j
                        pos = k
                        return [alpha, amplitude, pos]
                    else:
                        startingDir += 1

    # Feb14 High res moving scripts
    alphaL = [1e-3, 1e-4]
    amplitudeL = [3, 10]
    posL = [10, 30, 100, 300]
    startingDir = 202
    for i in alphaL:
        for j in amplitudeL:
            for k in posL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    pos = k
                    return [alpha, amplitude, pos]
                else:
                    startingDir += 1

    # Feb4 High Res stationary scripts
    alphaL = [1e-3, 1e-4]
    amplitudeL = [3, 10, 30]
    posL = [30, 60, 90]
    startingDir = 184
    for i in alphaL:
        for j in amplitudeL:
            for k in posL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    pos = k
                    return [alpha, amplitude, pos]
                else:
                    startingDir += 1

    # Jan 23 High Res stationary scripts
    alphaL = [1e-3, 1e-4]
    amplitudeL = [3, 10, 30]
    posL = [30, 60, 90]
    startingDir = 150
    for i in alphaL:
        for j in amplitudeL:
            for k in posL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    pos = k
                    return [alpha, amplitude, pos]
                else:
                    startingDir += 1

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

    if startingDir == 90:
        return [0.0001, 30, 30]

    # Dec 23 stationary scripts
    alphaL = [1e-3, 1e-4]
    amplitudeL = [3, 30]
    velocityL = [10, 30, 100, 300]
    startingDir = 130
    for i in alphaL:
        for j in amplitudeL:
            for k in velocityL:
                if startingDir == outputDir:
                    alpha = i
                    amplitude = j
                    vel = k
                    return [alpha, amplitude, vel]
                else:
                    startingDir += 1
    return [0,0,0]
