# _____ _
# | __ \ | |
# | |__) |__ _ __ _ __ ___ _ _| |_ ___ _ __ _ _
# | ___/ _ \ '__| '_ ` _ \| | | | __/ _ \ | '_ \| | | |
# | | | __/ | | | | | | | |_| | || __/ _ | |_) | |_| |
# |_| \___|_| |_| |_| |_|\__,_|\__\___| (_) | .__/ \__, |
# | | __/ |
# |_| |___/
## Written by John Lettman
## TCHS Computer Information Systems
## Professor Fuchs
#

# Scripting flags (alters permute())
OUTPUT_VERBOSE = True # Give additional state information.

ODOMETER_USEMINIFIED = False # Base-array on True, Class-array on False.
ODOMETER_INTEGERS = False # Integer counter. Turns off odometer mode.
ODOMETER_VISIBLE = False # Should we even use an odometer?

BENCHMARKING = True # Show timing.
COMPUTERSCIENCE = True # Start at 0, etc.
SYMBOLICNUMBERS = False # 2!, 3!, 4!, etc.

STRING_REFORMAT = True # All uppercased.

UNIQUE_SETINPUT = False # Unique chars input.
#UNIQUE_SETOUTPUT = False # Unique elements in the output list.


global noIncludes
# [!] Do not edit noIncludes.
noIncludes = False


#################
# Benchmarking. #
#################

if(BENCHMARKING):
    import time
    t0 = time.time()


###############################
# Import the proper odometer. #
###############################

if(not ODOMETER_INTEGERS):
    # Using Base-N odometer types.
    if(not ODOMETER_USEMINIFIED):
        # Use scientific source.
        try: import odometer as odometer
        except: noIncludes = True # Failed.


#########################
# Permutation function. #
#########################

def permute(a):
    """ Permutes all unique combinations of "a" as input."""
    """ Makes use of the yield keyword to render output before during execution."""

    ## Analyze user input differently for different types.
    ## Typically the function takes a list, but we can fix other inputs.
    if(isinstance(a, str)): # It's a string.
        if(STRING_REFORMAT):
            a = list(str(a).upper()) # str -> list
        else:
            a = list(a)
    elif(isinstance(a, int)): # It's an integer.
        a = list(str(a)) # int -> list(digits)
    elif(isinstance(a, list)):
       if(isinstance(a[0], list) and OUTPUT_VERBOSE):
           ## A list of lists...
           print("You're insane.")

    ## Make all elements unique.
    ## This is a rather "pythonic" way of processing such.
    if(UNIQUE_SETINPUT):
       set = {}
       map(set.__setitem__, a, [])
       a = set.keys()

    # a.sort() # Sort.

    ## Output the first input sorted.
    if(noIncludes or not ODOMETER_VISIBLE):
        yield list(a)
    elif(ODOMETER_INTEGERS):
        if(COMPUTERSCIENCE): intOdo = 0
        else: intOdo = 1
        yield list([intOdo, a])
    elif(not ODOMETER_USEMINIFIED):
        odo = odometer.Odometer(len(a))
        yield list([odo.GetValueSplit(" "), a])

    if(len(a) <= 1):
        ## No elements in "a" to permute!
        return # raise n00bProgrammerException

    i = 0 # Initialize our iterator var.
    first = 0 # We always start with 0 in an array!
    alen = len(a) # The length of our input.

    ## "alen" could also be used for the reference to the last element.

    while(True):
        i = alen - 1

        while(True):
            i -= 1 # i--

            if(a[i] < a[(i + 1)]):
                j = alen - 1

                while(not (a[i] < a[j])): j -= 1 # j--

                a[i], a[j] = a[j], a[i] # swap(a[j], a[i])
                t = a[(i + 1):alen]
                t.reverse()
                a[(i + 1):alen] = t

                # Output current.
                if(noIncludes or not ODOMETER_VISIBLE):
                    yield list(a)
                elif(ODOMETER_INTEGERS):
                    intOdo += 1
                    yield list([intOdo, a])
                elif(not ODOMETER_USEMINIFIED):
                    odo.Increment()
                    yield list([odo.GetValueSplit(" "), a])

                break # next.

            if(i == first):
                # Array of two, 2! = 2.
                # Best solution: reverse, add to orig, done.
                a.reverse()

                if(noIncludes or not ODOMETER_VISIBLE):
                    yield list(a)
                elif(ODOMETER_INTEGERS):
                    intOdo += 1
                    yield list([intOdo, a])
                elif(not ODOMETER_USEMINIFIED):
                    odo.Increment()
                    yield list([odo.GetValueSplit(" "), a])

                # End, we technically can't further permute this list.
                return


###########################
# PermuteToList function. #
###########################

def permuteToList(a):
    """ Accesses permute() and converts the yields to a list of lists."""
    product = [ ]
    for array in permute(a):
        # Constantly append to return until completion.
        product.append(array)

    return product


##############################
# PermuteFromInput function. #
##############################

# [!] permuteFromInput is non-functional currently.

##def permuteFromInput():
## """Accesses permute() with user data."""
## inputError = True
##
## while(inputError):
## inp = input("Enter a comma-separated list to permute:\n")
## if(inp == ""): inputError = True
## else: inputError = False
##
## return permute(str(inp).split(","))


######################
## Includes status. ##
######################

if(noIncludes and OUTPUT_VERBOSE):
    print("[!] Includes failed.\n\n")


#################
# Benchmarking. #
#################

if(BENCHMARKING):
    print("\nDone in: " + str((time.time() - t0)) + " seconds.")
