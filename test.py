

a = list(map(int, input().split()))
k=0
temp=[]
while k < len(a):
    p=0
    for l in range(len(a)):
        if a[k] > a[l]:
            p+=1
    k+=1
    temp.append(str(p))

print(" ".join(temp))


# 2

# def count(x, n, minimum):
#     if x == 0:
#         return 1
#     if x < minimum ** n:
#         return 0
#     return count(x - minimum ** n, n, minimum + 1) + count(x, n, minimum + 1)

# x, n = input().split()
# x, n = int(x), int(n)
# print(count(x, n, 1))