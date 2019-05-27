
.randomPowerLaw <-function(n, alpha=2, xmin=1) {

    # make required number of random input values between 0 and 1
    values = runif(n, 0, 1)

    # calculate power law function for these values and output results
    result = xmin * ((1-values) ^ (-1/(alpha-1)))
    return(result)
}

setMethod("randomPowerLaw",
          signature=signature(n="numeric"),
          .randomPowerLaw)
