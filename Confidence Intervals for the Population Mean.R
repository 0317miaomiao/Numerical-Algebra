set.seed(1)
sampledata<- rnorm(100,10,10)
#查看t.test中包含哪些列
ls(t.test(sampledata))

t.test(sampledata)$"conf.int"