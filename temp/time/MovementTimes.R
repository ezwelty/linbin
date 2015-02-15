# load data
d = read.csv('data/time/MovementTimes.csv', stringsAsFactors = F)

# format data
d = d[d$region != 'first read',]
d = d[order(d$tag),]
d$dtime = unclass(as.POSIXct(strptime(d$timestamp, format = "%m/%d/%y %H:%M")))
d$tag = as.character(d$tag)

# cut data at inconsistent antenna readings
tags = unique(d$tag)
for (tag in tags) {
  ind = d$tag == tag
  reads = nrow(d[ind,])
  if (reads > 1) {
    # check consistency
    r0 = substr(d[ind,'region'][1:(reads-1)], 2, 2)
    r1 = substr(d[ind,'region'][2:reads], 1, 1)
    errors = r0 != r1
    if (any(errors)) {
      d[ind, 'tag'] = paste(d[ind, 'tag'], as.character(cumsum(c(0, errors))), sep = '')
    }
  }
}

# process into interval data
n = data.frame(tag = character(), from = numeric(), to = numeric(), region = numeric(), stringsAsFactors = F)
tags = unique(d$tag)
# process into intervals
for (tag in tags) {
  temp = d[d$tag == tag,]
  reads = nrow(temp)
  if (reads > 1) {
    # parse into time intervals
    r0 = NULL
    for (i in 1:reads) {
      N = nrow(n)
      r1 = as.numeric(substr(temp$region[i],1,1))
      r2 = as.numeric(substr(temp$region[i],2,2))
      # first read
      if (is.null(r0)) {
        n[N+1,] = NA
        if (r1 == r2) {
          n$region[N+1] = r1
        } else {
          n$region[N+1] = r2
        }
        n$from[N+1] = temp$dtime[i]
        n$tag[N+1] = tag
        r0 = r2
      # last read
      } else if (i == reads) {
        n$to[N] = temp$dtime[i]
      # intermediate read
      } else {
        if (r1 == r2) {
          next
        } else {
          n$to[N] = temp$dtime[i]
          n[N+1,] = NA
          n$from[N+1] = temp$dtime[i]
          n$region[N+1] = r2
          n$tag[N+1] = tag
          r0 = r2
        }
      }
    } 
  } else {
    next
  }
}
n$ftag = substr(n$tag, 0, 15)

n_temp = n

# keep only fish with at least one region transition
n = n_temp
f = as.data.frame(table(n$ftag))
fish = as.character(f$Var1[f$Freq > 2])
ind = n$ftag %in% fish
n = n[ind,]

## keep only fish with no periods > 2 days
# tags = unique(n$tag)
# for (tag in tags) {
#   ind = n$tag == tag
#   temp = n[ind,]
#   if (any(temp$to - temp$from > 60*60*24*2)) {
#     n = n[!ind,]
#   }
# }

# FORMAT DATA
nd = import.data(
	file = n,
	fromField = "from",
	toField = "to",
	allow.nolength = FALSE
	)

# BIN DATA
binSizes = 60*60*(24/24)
binData = make.bins(
	rawFile = nd,
	binSizes = binSizes,
	factors = "region",
	factor.length = FALSE,
	equal.weight = TRUE,
	binRange = c(1217313360,1219178160)
	)

###

# bin fish by starting region
nr = n
ftags = unique(nr$ftag)
for (ftag in ftags) {
  ind = nr$ftag == ftag
  start_region = nr[ind, 'region'][1]
  tags = unique(nr[ind, 'tag'])
  for (tag in tags) {
    ind = nr$tag == tag
    temp = nr[ind,]
    temp$to[1] = temp$to[nrow(temp)]
    temp$region[1] = start_region
    nr = rbind(temp[1,], nr[!ind,])
  }
}
nr = import.data(
	file = nr,
	fromField = "from",
	toField = "to",
	allow.nolength = FALSE
	)
rbinData = make.bins(
	rawFile = nr,
	binSizes = binSizes,
	factors = "region",
	factor.length = FALSE,
	equal.weight = TRUE,
	binRange = c(1217313360,1219178160)
	)

#binData$r1n = binData$region.1 / (binData$region.1+binData$region.2)	
binData$r1n = binData$region.1 / (rbinData$region.1)

# PLOT DATA
barplot.grid(
	binFile = binData,
	fieldGroups = c("r1n"),
	binSizes = binSizes,
	dim = c(1,1),
	xLabels = 0
	)
# add daily vertical lines and time axis
xstart = unclass(as.POSIXct("2008-07-29"))[1]
xstop = unclass(as.POSIXct("2008-08-19"))[1]
days = seq(xstart, xstop, 3600*24)
for (d in days) {
  abline(v = d, col = 'grey')
}
axis(1, at = days[c(1,length(days))], labels = format(as.POSIXct(days[c(1,length(days))], origin = '1970-01-01'), '%b-%d'))

# 1217313360,1219178160
# 1217318400,1219178160
# n
ind = n$from >= 1217318400 & n$to <= 1219172700
motion = n[ind, c("ftag", "region", "from", "to")]
names(motion) = c("fish.id", "region", "from", "to")
motion = motion[order(motion$from, motion$to), ]
# nr
origin = nr[c("ftag", "region", "from", "to")]
names(origin) = c('fish.id', 'region', 'from', 'to')
origin = origin[order(origin$from, origin$to), ]
#
d = list(motion, origin)
fish = d
save('fish', file = 'data/fishmotion.rda')
