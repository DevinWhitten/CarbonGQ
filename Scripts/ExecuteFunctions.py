def Process_Arguments(array):
    #Precondition: Accepts the array of system arguments passed to execute.
    #Postcondition: Handles three conditions, no arguments, right arguments,
    #               and too many arguments. Just open sys.argv correctly

    if len(array) == 1: # No arguments, just return empty string.
        print "No keyword options. - Okay"
        return ""

    elif len(array) == 2: # Perfect, just grab the arguments
        return array[1]

    else:  # Too many arguments!
        print "Arguments not understood... ex. -v "
        print "Running in default mode ..."
        return ""
