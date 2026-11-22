test_that("duplicate measurement keys preserve exact times and input order", {
  id <- c("z","a","a","a","z","z","a","z","a")
  time <- c(1,0,-0,1,1,1+.Machine$double.eps,1+.Machine$double.eps,1,0)
  for (rows in list(integer(),1L,seq_along(id),rev(seq_along(id)),
                    c(5,1,4,2,3,7,9,6,8))) {
    key <- data.frame(id=id[rows],time=time[rows])
    expect_identical(.gpc_repeated_measurement(key$id,key$time),duplicated(key))
  }
})

test_that("interleaved unequal integer histories flag the later record", {
  rows <- data.frame(id=c("b","c","a","b","b","a","b"),time=c(0L,0L,0L,2L,1L,0L,2L))
  expect_identical(.gpc_repeated_measurement(rows$id, rows$time), duplicated(rows))
})
