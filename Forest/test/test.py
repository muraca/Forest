from subprocess import Popen,PIPE
import re
import time

tests = [2,4]
dimensions = [100, 500, 1000, 2500]
for i in tests:
    print("[+]TESTING with {} processes".format(i))
    for j in dimensions:
        print("DIM {}".format(j))
        entry = {}
        start = time.time()
        process = Popen(["./test.sh",str(i),str(j)], stdout=PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        end = time.time()
        t = end - start
        print("[{}] PROCESS TIME {}s".format(i,t))
