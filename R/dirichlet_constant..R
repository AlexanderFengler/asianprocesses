x <- rdirichlet(20, c(1,1,1) )

ddirichlet(x, c(1,1,1) )

a = 25
k = 3
x  = c(1 + (a/k),2 + (a/k), 1 + (a/k))
y  = c(a/k,a/k,a/k)
g_a = gamma(a)
g_n_a = gamma(4+a)

constant_griffiths = (prod(gamma(x))/prod(gamma(y)))*(prod(g_a)/prod(g_n_a))

constant_me = diri_norm(x)/diri_norm(y)
