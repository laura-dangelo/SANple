library(SANple)

n1 = 40
n2 = 30
n3 = 50
n4 = 60
set.seed(123)
sdd = 0.7
y1 = c(rnorm(n1,-5,sdd), rnorm(n2,0,sdd), rnorm(n3,5,sdd))
cl1 = c(rnorm(n1,-5,0), rnorm(n2,0,0), rnorm(n3,5,0))

y2 = c(rnorm(n1,-5,sdd), rnorm(n2,0,sdd), rnorm(n3,5,sdd))
cl2 = c(rnorm(n1,-5,0), rnorm(n2,0,0), rnorm(n3,5,0))

y3 = c(rnorm(n4,0,sdd), rnorm(n2,5,sdd), rnorm(n3,10,sdd))
cl3 = c(rnorm(n4,0,0), rnorm(n2,5,0), rnorm(n3,10,0))

y4 = c(rnorm(n1,0,sdd), rnorm(n4,5,sdd), rnorm(n3,10,sdd))
cl4 = c(rnorm(n1,0,0), rnorm(n4,5,0), rnorm(n3,10,0))

y = c(y1,y2,y3,y4)
cl = c(cl1,cl2,cl3,cl4)
group = c(rep(1,length(y1)), rep(2,length(y2)), rep(3,length(y3)), rep(4,length(y4)))
plot(y, col=group)

run1 = SANple::sample_CAM(3000, y, group, nclus_start = 5)
run1
clus1 = estimate_clusters(run1)
clus1
plot(run1, estimated_clusters = clus1)



run2 = SANple::sample_fiSAN(3000, y, group, nclus_start = 5, eps_beta = 0.3)
run2
clus2 = estimate_clusters(run2)
clus2
plot(run2, estimated_clusters = clus2, type = "ecdf")

### questo non va benissimo :(
run3 = SANple::sample_fiSAN_sparsemix(5000, y, group, nclus_start = 5, beta = 0.1)
run3
clus3 = estimate_clusters(run3)
clus3
plot(run3, estimated_clusters = clus3)
traceplot(run3, params = c("mu"), trunc_plot = 3)


run4 = SANple::sample_fSAN(3000, y, group, nclus_start = 5, beta = 1, alpha = 1)
run4
clus4 = estimate_clusters(run4)
clus4
plot(run4, estimated_clusters = clus4)


run5 = SANple::sample_fSAN_sparsemix(3000, y, group, nclus_start = 5, beta = 0.1, alpha = 0.1)
run5
clus5 = estimate_clusters(run5)
clus5
plot(run5, estimated_clusters = clus5)
