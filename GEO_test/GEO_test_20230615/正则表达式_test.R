#正则表达式练习
# 中括号[]表示选项，代表内部数据任意选择，比如[ATCG]，表示A，T，C，G四个字符随意选择。
# 如果是数字也是一个意思[1356],表示有1，2，5，6这个四个选项
# 如果嫌麻烦，可用-连接起始代表范围，
# 比如[A-Z]，代表大写的26个字母，
# 比如[a-z]，代表小写的26个字母，
# 比如[0-9]，代表从0到9的10个数字，
# 还可以混写
# 比如[0-9a-zA-Z]，代表数字和字母
# 
# 大括号{}代表重复次数，比如[ATCG]{3}表示从A,T,C,G四个字符中选择3个，那么这时候会产生ATC,ATG,ACG,TCG4种组合，在这里是为了形成三联密码子。
# 如果是两个数，用逗号隔开，代表范围，[0-9]{4,10} 表示产生4位数到10位数都可以。
# 如果两个数中的第二个缺失，比如[0-9]{4,}，代表4位数及以上。
# 
# 圆括号()代表成组，跟我们平常理解的意思一样，代表分组。实际上他的用处很多，但是分组，是核心功能。

strings <- c(
  "apple", 
  "219 733 8965", 
  "329-293-8753", 
  "Work: 579-499-7527; Home: 543.355.3679"
)

pattern <- "([1-9][0-9]{2})[- .]([0-9]{3})[- .]([0-9]{4})"

str_view(strings, pattern)

#str_detect，从strings中检测能否和pattern匹配，返回的是逻辑值，有点像grepl
str_detect(strings, pattern)
# [1] FALSE  TRUE  TRUE  TRUE

#str_locate，从strings中检测能否和pattern匹配，返回匹配的起止位置
str_locate(strings,pattern)

#3.取回
#第一，str_subset返回的是匹配到的原始条目
str_subset(strings, pattern)

#第二，str_extract返回的是匹配到的模式
str_extract(strings, pattern)

#第三，str_match返回的是数据框
str_match(strings, pattern)
#第一列是str_extract的数据，后面依次是括号中的内容，模式中有多少个(),就返回多少列。本例中有三对小括号。

#str_replace从strings中能够精确匹配pattern的内容，并替换为''XXX-XXX-XXXX''(此处自定义)
str_replace(strings, pattern, "123-456-7894")



# 问号?表示0次或者是1次，因为这是一个生存或是毁灭的问题
# 加号+表示1次或者多次，把加号和1联系起来，用医(1)院来记忆
# 星号*表示0次或者多次，把星号和零联系起来，用零(0)星来记忆




