from subprocess import Popen,PIPE
import re

tests = [1,2,4]
dimensions = [100, 1000, 10000, 100000]
for i in tests:
    print("[+]TESTING with {} processes".format(i))
    for j in dimensions:
        print("DIM {}".format(j))
        entry = {}
        process = Popen(["./test.sh",str(i),str(j)], stdout=PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        output = output.splitlines()
        for x in output:
            temp = x.decode().split(' ')[-1]
            code = re.search(r"\[([0-9]+)\]",x.decode()).group(1)
            entry[int(code)] = int(temp)
        roottime = entry[0]/1000000
        print("[{}] PROCESS TIME {}s".format(i,roottime))
