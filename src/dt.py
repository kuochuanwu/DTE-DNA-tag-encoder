# -*- coding: utf8 -*-
import copy, math, os



class ReadFile:
    '''注意，判斷不同大小序列跳出error還沒寫'''
    path = os.path.split(os.path.realpath(__file__))[0] + "//data//test.fas"
    # path = os.path.split(os.path.realpath(__file__))[0] + "//data//Brassicaceae_cuted.fas"
    
    comment = []
    sequence = []
    sampleSize = 0
    sequenceLength = 0
    try:
        f = open(path, "r")
        for line in f:
            #getting comment of data
            if line[0] is ">":
                comment.append(line.upper().strip("\n"))
            else:
                sequence.append(line.upper().strip("\n"))
            #getting attribute of data
        sampleSize = len(sequence)
        sequenceLength = len(sequence[0])  
    finally:
        f.close()
     

class TreeStructure:
    '''
        利用節點陣列表示法(鏈結串列表示法 linked-list representation)：
        left | data | right
            
              A(0)
           /       \
          B(1)     C(2)
           \      /
            D(3) E(4)
            
        |   index   |   left    |   data    |   right
            0           1           A           2
            1           -1          B           3
            2           4           C           -1
            3           -1          D           -1
            4           -1          E           -1
    '''
    
    '''|Index|Left|Set(index)|sequence(reduced)|node|SNP|Right|'''
    #放樣本集合
    set_index = []
    index = 0
    leftChild = 0
    sequence = []
    node = 0
    rightChild = 0




'''初始類別'''
#讀檔
DATA = ReadFile
#TREE = TreeStructure

'''初始參數'''
data = ["" for i in range(len(DATA.sequence))]
data_index = []

#放樹狀結構
'''|Index|Left|Set(index)|node|Right|'''
forest = []
INDEX = 0
INDEXLEAF = 0
indexLeaf = 0


'''資料簡化(去掉一樣的與忽略N)'''
def formatData():
    tempChar = ''
    for i in range(len(DATA.sequence[0])):
        tempChar = DATA.sequence[0][i]
        same = True
        for j in range(len(DATA.sequence)):
            if tempChar != DATA.sequence[j][i]:
                N_appear = False
                for j in range(len(DATA.sequence)):
                    if 'N' == DATA.sequence[j][i]:
                        N_appear = True
                        break
                if N_appear == False:
                    same = False
                break
        if same is False:
            data_index.append(i)
    for i in range(len(data_index)):
        for j in range(len(DATA.sequence)):
            data[j] += DATA.sequence[j][data_index[i]]


'''找最接近目標值(樣本中間數)'''
def closeTarget(num, target):
    mid = 9999999
    '''
        result = [
            接近目標值的核苷酸索引值 0:A, 1:C, 2:G, 3:T
            接近目標值的差,
            核苷酸分佈數量區間,
            
        ]
    '''
    
    
    result = [0, 0, 0]
    result_index = 0
    for i in range(len(num)):
        temp = math.fabs(target - num[i])
        if mid > temp:
            mid = temp
            result[0] = i
            result[1] = temp
        if num[i] != 0:
            result[2] += 1
        
    if result[2] == 1:
        result[2] = 0
    elif result[2] == 2:
        result[2] = 1.0
    elif result[2] == 3:
        result[2] = 0.66
    elif result[2] == 4:
        result[2] = 0.33
    
    # print num, result[0], "\t", 
    return result



def toSeparate(tree):
    global INDEX
    global indexLeaf
    ACGT = []
       
    for i in range(len(data[0])):
        acgt = [0,0,0,0]
        for j in range(len(tree[2])):
            if data[tree[2][j]][i] is "A":
                acgt[0] += 1
            elif data[tree[2][j]][i] is "C":
                acgt[1] += 1
            elif data[tree[2][j]][i] is "G":
                acgt[2] += 1
            elif data[tree[2][j]][i] is "T":
                acgt[3] += 1
        ACGT.append(acgt)
    #最佳分群分數
    max_score = 0
    #最佳分群核苷酸位點
    max_index = 0
    #分群位點之核苷酸
    max_acgt = ''
    score = 0.0
    temp_acgt = 0
    mid = float(len(tree[2]))/2.0
    for i in range(len(ACGT)):
        result = closeTarget(ACGT[i], mid)
        score = (mid - result[1]) / mid + result[2]
        # print score
        if max_score < score:
            max_score = score
            max_index = i
            temp_acgt = result[0]
    # print         
        
    if temp_acgt == 0:
        max_acgt = 'A'
    elif temp_acgt == 1:
        max_acgt = 'C'
    elif temp_acgt == 2:
        max_acgt = 'G'
    elif temp_acgt == 3:
        max_acgt = 'T'
    

    forest[INDEX][3] = max_index
    
    #print forest
    #print tree.set_index
    #print "forest[INDEX].leftChild: ", forest[INDEX][1]
    #print "forest[INDEX].rightChild: ", forest[INDEX][4]
    
    '''若一節點超過2個集合，向下延伸兩個節點'''
    forest.append([0,0,[],0,0])
    forest.append([0,0,[],0,0])
        
    '''放集合'''
    if len(tree[2]) == 2:
        '''若一節點剩2個集合，向下延伸左右節點'''
        forest[forest[INDEX][1]][2].append(tree[2][0])
        forest[forest[INDEX][4]][2].append(tree[2][1])
    else:
        left = 0
        right = 0
        '''以最接近中位數目標做為分兩類依據'''
        for i in range(len(tree[2])):
            #print "data[tree.set_index[i]][max_index]: ", data[tree.set_index[i]][max_index]
            #print "max_acgt: ", max_acgt
            if data[tree[2][i]][max_index] is max_acgt:
                #print "TEST"
                '''Left'''
                forest[forest[INDEX][1]][2].append(tree[2][i])
                left += 1
            else:
                '''Right'''
                forest[forest[INDEX][4]][2].append(tree[2][i])
                right += 1
        #if left == 0 or right == 0:
        #    forest[forest[INDEX][1]][2] = []
        #    forest[forest[INDEX][4]][2] = []
        #    forest[forest[INDEX][1]][2].append(tree[2][0])
        #    for i in range(len(tree[2])):
        #        forest[forest[INDEX][4]][2].append(tree[2][i])
    
    # print "INDEX:", INDEX
    # for i in range(len(forest)):
    #     print forest[i]
    # print
    # for i in range(len(forest)):
    #     print forest[i]
    
    #print forest[forest[INDEX].leftChild].leftChild
    #print "forest[INDEX].leftChild: ", forest[INDEX].leftChild
    
    '''設定下個葉子的指標'''
    #左葉
    if len(forest[forest[INDEX][1]][2]) == 1:
        forest[forest[INDEX][1]][1] = -1
        forest[forest[INDEX][1]][4] = -1
    elif len(forest[forest[INDEX][1]][2]) > 1:
        forest[forest[INDEX][1]][1] = indexLeaf + 1
        forest[forest[INDEX][1]][4] = indexLeaf + 2
        indexLeaf += 2
    #右葉
    if len(forest[forest[INDEX][4]][2]) == 1:
        forest[forest[INDEX][4]][1] = -1
        forest[forest[INDEX][4]][4]= -1
    elif len(forest[forest[INDEX][4]][2]) > 1:
        forest[forest[INDEX][4]][1] = indexLeaf + 1
        forest[forest[INDEX][4]][4] = indexLeaf + 2
        indexLeaf += 2


formatData()

print("Sample size: " + str(DATA.sampleSize))
print("Length of sequences: " + str(DATA.sequenceLength))
print("After format length of sequences: " + str(len(data_index)))

print("DATA:")
for i in data:
    print (i)
print()
print("short data index:")
print(str(data_index))
print()

forest.append([0,0,[],0,0])
forest[INDEX][1] = 1
forest[INDEX][4] = 2
indexLeaf = indexLeaf + 2;
for i in range(len(data)):
   forest[INDEX][2].append(i)


'''資料存取'''
while(True):
    forest[INDEX][0] = INDEX
    #print INDEX, ":", len(forest[INDEX][2])
    #print "indexLeaf:", indexLeaf
    #print forest[INDEX][2]
    
    #print "INDEX: ", INDEX, "\tlen(forest): ", len(forest), "\tindexLeaf", indexLeaf
    #print forest[INDEX]
    if len(forest[INDEX][2]) > 1:
        toSeparate(forest[INDEX])
    if INDEX == indexLeaf:
        break
    if INDEX > 2**DATA.sampleSize:
        break
    INDEX += 1

#顯示樹狀結構
for i in range(len(forest)):
    print (str(forest[i]))

node = []
barcode = []
for i in range(len(forest)):
    if len(forest[i][2]) > 1:
        repeat = False
        for j in range(len(node)):
            if(forest[i][3] == node[j]):
                repeat = True
        if not repeat:
            node.append(forest[i][3])
node.sort()


for i in range(len(data)):
    str = ""
    for j in node:
        str += data[i][j]
    barcode.append(str)

print("")
print("data_index")
print(data_index)
for j in node:
    print(data_index[j], end=" ")
print("")
print("")
    
# for i in range(len(barcode)):
#     print "%d:\t" %i, barcode[i]

# for i in range(len(barcode)):
#     for j in barcode[i]:
#         print "%s\t" % j,
#     print("")
# print("")

for i in range(len(barcode)):
    print(DATA.comment[i])
    print(barcode[i])

print(len(barcode[0]))

