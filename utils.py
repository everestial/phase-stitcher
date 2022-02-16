import resource


def print_authorship():
    print("#######################################################################")
    print("        Welcome to phase-Stitcher version %s       " % 1.2)
    print("  Author: kiran N' bishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) ")
    print("#######################################################################")
    print()


# ''' to monitor memory usage. '''
def current_mem_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
