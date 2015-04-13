__author__ = 'Dmitrii'


import sys
from collections import defaultdict

file_first_solution = open(str(sys.argv[1]), 'r')
file_second_solution  = open(sys.argv[2], "r")


first_solution = {}
second_solution = {}


for row in file_first_solution.readlines():
    splitted_row = row.split()
    first_solution[splitted_row[0]] = int(splitted_row[1])

for row in file_second_solution.readlines():
    splitted_row = row.split()
    second_solution[splitted_row[0]] = int(splitted_row[1])


blast_file  = open(sys.argv[3], "r")




lengthes = {}
lengthes_just = {}

lengthes_short = []
lengthes_medium = []
lengthes_long = []


answer = defaultdict(int)

f3 = open("summary.txt", 'w')
for row in blast_file.readlines():
    splitted_row = row.split()
    if splitted_row[0] not in lengthes and splitted_row[0] not in lengthes_just and float(splitted_row[2]) > 95.0 and splitted_row[3] == splitted_row[7] and splitted_row[6] == "1" and splitted_row[4] == "0" and splitted_row[5] == "0" :
        lengthes[splitted_row[0]] = splitted_row[3]
        if int(splitted_row[7]) < 250:
            lengthes_short.append(splitted_row[0])
        if int(splitted_row[7]) >= 250 and int(splitted_row[7]) < 1000 :
            lengthes_medium.append(splitted_row[0])
        if int(splitted_row[7]) >= 1000 :
            lengthes_long.append(splitted_row[0])
    if splitted_row[0] in lengthes and float(lengthes[splitted_row[0]]) ==  float(splitted_row[3]) and float(splitted_row[2]) > 95.0:
        answer[splitted_row[0]] += 1
    lengthes_just[splitted_row[0]] = 1
f1 = open("difference_magnolia.txt", 'w')

difference_first = 0
difference_first_short = 0
difference_first_medium = 0
difference_first_long = 0
true_first_short = 0
true_first_medium = 0
true_first_long = 0

for i in answer.keys():
    if int(float(lengthes[i])) > 0:
        if i in lengthes_short and i in answer and i in first_solution:
            difference_first_short += abs(answer[i] - first_solution[i])
            if abs(answer[i] - first_solution[i]) == 0:
                true_first_short += 1
        if i in lengthes_medium and i in answer and i in first_solution:
            difference_first_medium+= abs(answer[i] - first_solution[i])
            if abs(answer[i] - first_solution[i]) == 0:
                true_first_medium += 1

        if i in lengthes_long and i in answer and i in first_solution:
            difference_first_long += abs(answer[i] - first_solution[i])
            if abs(answer[i] - first_solution[i]) == 0:
                true_first_long += 1

        if i not in first_solution:
            difference_first += answer[i]
            f1.write(i + " " + str(answer[i]) + "\n")
        else:
            difference_first += abs(answer[i] - first_solution[i])
            f1.write(i + " " + str(answer[i] - first_solution[i]) + "\n")

f3.write("Summary for Magnolia")
f3.write("Total difference " + str(difference_first) + "\n")
f3.write("Total difference short <250bp " + str(difference_first_short) + "\n")
f3.write("True short <250bp " + str(true_first_short) + "\n")


f3.write("Total difference medium 250-1000bp " + str(difference_first_medium) + "\n")
f3.write("True medium 250-1000bp " + str(true_first_medium) + "\n")

f3.write("Total difference long >1000 bp " + str(difference_first_long) + "\n")
f3.write("True long >1000bp " + str(true_first_long) + "\n")


f2 = open("difference_CNVera.txt", 'w')
f3.write("=============================")

difference_second = 0
difference_second_short = 0
difference_second_medium = 0
difference_second_long = 0
true_second_short = 0
true_second_medium = 0
true_second_long = 0

for i in answer.keys():
    if int(float(lengthes[i])) > 0:
        if i in lengthes_short and i in answer and i in first_solution:
            difference_second_short += abs(answer[i] - second_solution[i])
            if abs(answer[i] - second_solution[i]) == 0:
                true_second_short += 1
        if i in lengthes_medium and i in answer and i in first_solution:
            difference_second_medium+= abs(answer[i] - second_solution[i])
            if abs(answer[i] - second_solution[i]) == 0:
                true_second_medium += 1

        if i in lengthes_long and i in answer and i in first_solution:
            difference_second_long += abs(answer[i] - second_solution[i])
            if abs(answer[i] - second_solution[i]) == 0:
                true_second_long += 1

        if i not in second_solution:
            difference_second += answer[i]
            f2.write(i + " " + str(answer[i]) + "\n")
        else:
            difference_second += abs(answer[i] - second_solution[i])
            f2.write(i + " " + str(answer[i] - second_solution[i]) + "\n")

f3.write("Summary for CNVera")
f2.write("Total difference " + str(difference_second) + "\n")
f2.write("Total difference short <250bp " + str(difference_second_short) + "\n")
f2.write("True short <250bp " + str(true_second_short) + "\n")

f2.write("Total difference medium 250-1000bp " + str(difference_second_medium) + "\n")
f2.write("True medium 250-1000bp " + str(true_second_medium) + "\n")

f2.write("Total difference long >1000 bp " + str(difference_second_long) + "\n")
f2.write("True long >1000bp " + str(true_second_long) + "\n")
f3.write("==============================")

f3.write("Total summary")
f2.write("Total short - " + str(len(lengthes_short)) + "\n")
f2.write("Total medium - " + str(len(lengthes_medium)) + "\n")
f2.write("Total long -  " + str(len(lengthes_long)) + "\n")

