sfile      <- "sfile"
kernel_out <- "kernel"
rm_out     <- "rm_reactions"
rmlist_out <- "rm_list"

# --- DON'T TOUCH ---

library('pracma')
options(digits = 13)
m <- as.matrix(read.table(sfile))
k <- nullspace(m)
reactions <- nrow(k)
rm('m')
temp <- t(apply(k, 1, abs))
ix <- apply(temp, 1, which.max)
vals <- cbind(1:nrow(k), ix)
testk <- k/k[vals]
testk[testk == 0] <- NA
rm(list = c('temp', 'ix', 'vals'))
rmlist <- rep(0, nrow(testk))
removed <- c()
for (i in 1:(nrow(testk) - 1)) {
    print(i)
  if (!(i %in% removed)) {
    for (j in (i + 1):nrow(testk)) {
      if (!(j %in% removed)) {
        if (isTRUE(all.equal(testk[i,], testk[j,]))) {
          rmlist[j] <- i
          removed <- c(removed, j)
        }
      }
    }
  }
}
removed <- unique(sort(removed))
rm(list = c('testk', 'i', 'j'))
kernel <- k[-removed,]
rm('k')
removed <- removed - 1
rmlist <- rmlist - 1
write.table(kernel, file=kernel_out, row.names = F, col.names = F, sep=" ")
write.table(removed, file=rm_out, row.names = F, col.names = F)
write.table(rmlist, file=rmlist_out, row.names = F, col.names = F)
