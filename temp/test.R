### Simple Stupid

from = c(0,5,20,30,45)
to = c(10,15,30,40,50)
x = c(1,2,3,NA,1)
f = c('a','b',NA,'b','a')
e = event.table(from, to, x, f)
b = seq.bin.event.table(event.coverage(e), 4)
eb = sample.events(e, b, list(sum, 'x', na.rm = T), list(sum, 'x', by = 'f', na.rm = T))
plot.events(eb, 3)


df <- data.frame(grp = c("a", "a", "b", "b", "c", "c", "c"),
                 val = rnorm(7))
df %>% group_by(grp) %>% summarise(total = sum(val))

###
## Meta functions
citation("linbin")
packageDescription("linbin")
packageVersion("linbin")

## Initialization functions
verbose = TRUE
file = "../test.csv"
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M', verbose = verbose)

# coerce and validate (fail)
EventTable(c(1,2), data.frame(c(1),c(2)), verbose = verbose)  # multiple columns
EventTable(data.frame(ASD = c(1,2,4,5)), data.frame(c(1,2)), verbose = verbose)  # coerced to same length
EventTable(FROM = numeric(), TO = numeric(), verbose = verbose)	# Empty
EventTable(FROM=c('a','b'), TO=1, verbose = verbose)	# Not numeric
EventTable(FROM=1, TO=NA, verbose = verbose)	# Not finite (length 1) -- throws error as non numeric (NA alone is logical)
EventTable(FROM=c(1,2), TO=c(10,NA), verbose = verbose)	# Not finite (length > 1)
EventTable(FROM=c(1,2), TO=c(10,1), data.frame(FROM=c(4,5)), verbose = verbose)	# Repeating columns
is.EventTable(c(1), verbose = verbose)	# wrong class
is.EventTable(data.frame(), verbose = verbose)  # wrong class
import.EventTable('notafile', verbose = verbose)	# incorrect filename
import.EventTable(1, verbose = verbose)	# filename not character
import.EventTable(file, sep = ",", verbose = verbose)	# missing FROM, TO
import.EventTable(file, sep = ",", fromCol = 'FACTOR', toCol = 'TO_M', verbose = verbose)	# field not numeric
import.EventTable(file, sep = ",", fromCol = 'FROM', toCol = 'TO_M', verbose = verbose)	# field not found
as.EventTable(data.frame(FROM=1,TO=2,FROM_M=3),fromCol="FROM_M", verbose = verbose) # FROM field already exists

# coerce and validate (success)
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M', verbose = verbose)
is.EventTable(events, verbose = verbose)
validate.EventTable(events, validate.sort = TRUE, verbose = verbose)

# validate bins
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')
validate.binEventTable(events, verbose = verbose) # Not bins (points)
validate.binEventTable(EventTable(c(0,10,15),c(10,20,30)), verbose = verbose) # Not bins (overlapping)
validate.binEventTable(EventTable(c(0,10,20),c(10,20,30)), verbose = verbose) # bins (manual)


###
## Sorting

temp = data.frame(FROM=c(1,1,1,3,2), TO=c(1,4,10,5,20))	# not sorted on FROM
class(temp) = c("EventTable", class(temp))
validate.EventTable(temp, verbose = verbose)

temp = data.frame(FROM=c(1,1,1,2,3), TO=c(1,10,4,5,20))	# not sorted on TO
class(temp) = c("EventTable", class(temp))
validate.EventTable(temp, verbose = verbose)

temp = data.frame(FROM=c(1,1,1,2,20), TO=c(1,4,10,5,3))	# Switched FROM and TO, without affecting sort
class(temp) = c("EventTable", class(temp))
validate.EventTable(temp, verbose = verbose)
temp = as.EventTable(temp)  # Now sorted!
validate.EventTable(temp, verbose = verbose)

temp = data.frame(FROM=c(1,1,1,2,3), TO=c(1,4,10,5,20))  # Sorted!
class(temp) = c("EventTable", class(temp))
validate.EventTable(temp, verbose = verbose)

###
## Logical finds

# points
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')
pointpos = find.points(events, verbose = verbose)
events[pointpos,]

# lines
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')
linepos = find.lines(events, verbose = verbose)
events[linepos,]

# intersecting
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')
find.intersecting(c(1,2,3), events, verbose = verbose) # length > 2
find.intersecting(events, events, verbose = verbose) # nrow > 1
find.intersecting('a', events, verbose = verbose) # wrong class
find.intersecting(c(1,2), 'a', verbose = verbose) # wrong class
events[find.intersecting(c(200,300), events, closed.lines = TRUE, closed.points = TRUE, verbose = verbose),]
events[find.intersecting(c(200,300), events, closed.lines = FALSE, closed.points = TRUE, verbose = verbose),]
events[find.intersecting(c(200,300), events, closed.lines = FALSE, closed.points = FALSE, verbose = verbose),]

###
## EventTable Properties
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')

range.EventTable(events, verbose = verbose)
coverage.EventTable(events, closed = TRUE, verbose = verbose)
coverage.EventTable(events, closed = FALSE, verbose = verbose)
gaps.EventTable(events, verbose = verbose)
gaps.EventTable(EventTable(0,100), verbose = verbose) # no gaps

###
## EventTable Operations
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')

transform.EventTable(events, 'a', 100, verbose = verbose)  # wrong class
transform.EventTable(events, c(1,2), 100, verbose = verbose)  # wrong length
transform.EventTable(events, 2, 500, verbose = verbose)

cut.EventTable(events, 'a', totals = "COUNT", verbose = verbose) # wrong class
cut.EventTable(events, 69, totals = "COUNT", verbose = verbose) # single
cut.EventTable(events, c(0,100,200,300), totals = "COUNT", verbose = verbose) # multiple
cut.EventTable(events, EventTable(FROM=50,TO=100), totals = "COUNT", verbose = verbose) # event table

crop.EventTable(events, EventTable(FROM=c(0),TO=c(25)), totals = "COUNT", verbose = verbose)
crop.EventTable(events, EventTable(FROM=c(0,300),TO=c(100,400)), totals = "COUNT", verbose = verbose)

seq.binEventTable('a', bincount = 10, verbose = verbose) # wrong type
seq.binEventTable(1, bincount = 10, verbose = verbose) # scalar coverage
seq.binEventTable(c(1,1), bincount = 10, verbose = verbose) # zero range
seq.binEventTable(c(100,200), bincount = 10, verbose = verbose) # vector
seq.binEventTable(EventTable(100,200), bincount = 10, verbose = verbose) # events
seq.binEventTable(events, bincount = 10, verbose = verbose) # coverage
seq.binEventTable(range.EventTable(events), bincount = 10, verbose = verbose) # range
seq.binEventTable(c(0,100), bincount = 10, binlength = 1, verbose = verbose) # both defined
seq.binEventTable(c(0,100), bincount = 10, binlength = 1, verbose = verbose) # both defined
seq.binEventTable(c(0,100), verbose = verbose) # neiter defined
seq.binEventTable(c(0,100), binlength = 75, adaptive.binlength = T, verbose = verbose) # adaptive
gaps.EventTable(events)
gaps.EventTable(seq.binEventTable(events, binlength = 75, adaptive.binlength = T, verbose = verbose)) # adaptive

#################################
## Speed Tests
# data.table a little faster (especially on rbind such things), but data.matrix way faster!

# grow events
events = import.EventTable(file, sep = ",", fromCol = 'FROM_M', toCol = 'TO_M')
events2 = as.data.table(events)

new = events[rep(1:nrow(events), times=5000), ]
new2 = as.data.table(new)

system.time(for (i in 1:10) as.EventTable(new))
system.time(for (i in 1:10) as.EventTable(new2))

big = as.EventTable(new)
big2 = as.EventTable(new2)

system.time(for (i in 1:10) find.intersecting(c(0,450), big, closed.lines = TRUE, closed.points = TRUE, verbose = verbose))
system.time(for (i in 1:10) find.intersecting(c(0,450), big2, closed.lines = TRUE, closed.points = TRUE, verbose = verbose))

system.time(for (i in 1) coverage.EventTable(big))
system.time(for (i in 1) coverage.EventTable(big2))

system.time(for (i in 1:10) rbind(big,events))
system.time(for (i in 1:10) rbind(big2,events2))

big3 = as.matrix(big2)
events3 = as.matrix(events2)
system.time(for (i in 1:10) rbind(big3,events3))

N = 1
library(plyr)
system.time(for (i in 1:N) { 
  temp = as.data.frame(table(big$FROM,big$TO))
  temp = temp[temp$Freq > 0,]
})
system.time(for (i in 1:N) { 
  temp = ddply(big, .(FROM,TO), summarise, count = length(FROM))
})

#########################
## Complete (working) Test

## Meta functions
citation("linbin")
packageDescription("linbin")
packageVersion("linbin")

## Load and verify data
v = TRUE
# create
test = data.frame(from = c(1,2,3), to = c(2,3,4), data = c(0.1,0.2,0.3))
events = EventTable(FROM = test$from, TO = test$to, data = test$data, verbose = v)
events = as.EventTable(test, fromCol = "from", toCol = "to", verbose = v)
# import
events = import.EventTable('test.csv', fromCol = 'FROM_M', toCol = 'TO_M', sep = ',', verbose = v)
# verify
is.EventTable(events, verbose = v)
validate.EventTable(events, verbose = v)
validate.binEventTable(events, verbose = v)

## Logical indices
find.points(events, verbose = v)
find.lines(events, verbose = v)
find.intersecting(c(100,200), events, verbose = v)

## Measurements
range.EventTable(events, verbose = v)
coverage.EventTable(events, verbose = v)
gaps.EventTable(events, verbose = v)
#totlength.EventTable(events, verbose = v)
#ncoverage.EventTable(events, verbose = v)

## Manipulations
transform.EventTable(events, scale = 2, translate = 100, verbose = v)
cut.EventTable(events, c(0,100,200), totals = c('COUNT'), verbose = v)
crop.EventTable(events, as.EventTable(data.frame(FROM=0,TO=100)), totals = c('COUNT'), verbose = v)

## Building Bins
seq.binEventTable(range.EventTable(events), bincount = 10, verbose = v)
seq.binEventTable(range.EventTable(events), binlength = 50, verbose = v)
seq.binEventTable(range.EventTable(events), binlength = 50, adaptive.binlength = T, verbose = v)
gaps.EventTable(seq.binEventTable(coverage.EventTable(events), bincount = 10, verbose = v))
gaps.EventTable(seq.binEventTable(coverage.EventTable(events), binlength = 10, verbose = v))
gaps.EventTable(seq.binEventTable(coverage.EventTable(events), binlength = 10, adaptive.binlength = T, verbose = v))
validate.binEventTable(seq.binEventTable(coverage.EventTable(events), binlength = 10, adaptive.binlength = T, verbose = v))

## Binning
bins = build.binEventTable(c(0,500), binCount = 5, verbose = v)
overlaps = find.binEventOverlaps(events, bins)