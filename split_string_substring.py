#!usr/bin/python
# Hacker rank problem's solution. 
'''
Question:
-------------
String s is split into n/k = 9/3 = 3 equal parts of length k = 3. 
We convert each t[i] to u[i] by removing any subsequent occurrences non-distinct characters in t[i]:

1. t[i] = "AAB" --> u[i] = "AB"
2. t[i] = "CAA" --> u[i] = "CA"
3. t[i] = "ADA" --> u[i] = "AD"

We then print each  on a new line.
'''


def merge_the_tools(string, k):
    # your code goes here
    n = len(string)
    if(n % k != 0):
        exit
    else:
        seg_len = k
        res = []
        
        for i in range(0,len(string)-seg_len+1,seg_len):
            res.append(string[i:i+seg_len])

        for x in res:
            tmp = []
            res = []
            tmp = list(x)
            #print(''.join(set(tmp)))
            final = []
            [final.append(item) for item in tmp if item not in final]
            print(''.join(final))
        

if __name__ == '__main__':
    string, k = input(), int(input())
    merge_the_tools(string, k)
