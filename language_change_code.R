

# set global options

library(fmsb) # a library to compute Negelkerke R2 values


# function load and format the input datasets

attestations = function(file) {
    current.table = read.table(file)
    get.attestations = paste(current.table[,3], current.table[,4], sep="_")
    table(get.attestations)
}


# load the datastes

abo = attestations("data/abo.txt")
albo = attestations("data/albo.txt")
barzo = attestations("data/barzo.txt")
bardzo = attestations("data/bardzo.txt")
bych = attestations("data/bych.txt")
bym = attestations("data/bym.txt")
inny = attestations("data/inny.txt")
inszy = attestations("data/inszy.txt")
stopniowanie_na = attestations("data/stopniowanie_na.txt")
stopniowanie_naj = attestations("data/stopniowanie_naj.txt")
wiekszy = attestations("data/wiekszy.txt")
wietszy = attestations("data/wietszy.txt")
wszytko = attestations("data/wszytko.txt")
wszystko = attestations("data/wszystko.txt")
# the format for -ir-/-ier- is different:
ir = table(read.table("data/tmp_ir.txt"))
ier = table(read.table("data/tmp_ier.txt"))

# additional datasets
load("data/bardzo_ONLY.RData")
load("data/barzo_ONLY.RData")
load("data/bardziej_ONLY.RData")
load("data/barziej_ONLY.RData")
load("data/bychmy.RData")
load("data/bysmy.RData")


# custom colors for the plots
my.blue = rgb(0.15, 0.45, 0.96)
my.green = rgb(0.15, 0.85, 0.27, 0.7)
my.red = rgb(0.92, 0.3, 0.3, 0.6)
my.grey = rgb(0,0,0,.6)
my.orange = rgb(1,0.5,0.2,0.6) 
my.teal = rgb(0, 0.5, 0.5, 0.7) 
my.violet = rgb(0.75, 0.25, 0.82, 0.7)





# a function to compute the proportions of the old and the new form:
# its name comes from the Polish word 'uporządkuj' ('make order');
# the function takes the arguments 'zakres.chrono' (the subcorpus size)
# and 'krok' (lag, or the overlap for the moving window)

uporzadkuj = function(dawne, nowe, 
                      zakres.chrono = 20,
                      krok = 10,
                      zacznij = 1350, 
                      zakoncz = 1960) {

    dawne = dawne[sort(union(names(dawne), names(nowe)))]
    names(dawne) = sort(union(names(dawne), names(nowe)))
    dawne[is.na(dawne)] = 0

    nowe = nowe[sort(union(names(dawne), names(nowe)))]
    names(nowe) = sort(union(names(dawne), names(nowe)))
    nowe[is.na(nowe)] = 0
    
	results.all = c()
	shift.dates.all = c()
	skoki.chrono = seq(0, (zakres.chrono-krok), (zakres.chrono / (zakres.chrono/krok)) )
	
	for(i in skoki.chrono) {
			cezury = seq(zacznij, zakoncz, zakres.chrono) + i
			daty = as.numeric(gsub("^([0-9]{4}).+", "\\1", names(dawne)))
			przedzialy.czasowe = findInterval(daty, cezury)
			
			x = c()
			y = c()
			
			for(j in 1: length(cezury)) {
				x[j] = sum(nowe[przedzialy.czasowe == j])
				y[j] = sum(dawne[przedzialy.czasowe == j])
			}
			
			results = x / (x+y)
			
			# adjusting some values if dividing by O occurred
			results[is.nan(results)] = NA
			
			# shifting the scale so that the value between two dates is reached
			shift.dates = cezury + ((cezury[2] - cezury[1]) / 2)
			# getting rid of the last value
			#shift.dates = shift.dates[-length(shift.dates)]
			
			results.all = c(results.all, results)
			shift.dates.all = c(shift.dates.all, shift.dates)
	}
	
	# ordering the results
	results.all = results.all[order(shift.dates.all)]
	shift.dates.all = sort(shift.dates.all)
	names(results.all) = shift.dates.all
	
	# getting rid of NA values
	results.all = results.all[!is.na(results.all)]
return(results.all)
}









# Fig. 1, including model inference
# dataset here being the change więtszy > większy

#
piotrowski = uporzadkuj(wietszy, wiekszy, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim = c(0,1), xlab = "year", ylab = "proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
piotrowski1 = uporzadkuj(wietszy, wiekszy, zakres.chrono = 1, krok = 1)
model1 = glm(piotrowski1 ~ as.numeric(names(piotrowski1)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col = my.blue, lwd = 3)
lines(names(piotrowski1), model1$fitted, type="l", col = my.red, lwd = 2, lty = 2)
legend("topleft", c("więtszy > większy", "alternative model"), 
       text.col = c(my.blue, my.red),
       bty = "n",
       lty = c(1, 2), col = c(my.blue, my.red), lwd = 3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))








# Fig. 2. the change barzo > bardzo
#
piotrowski = uporzadkuj(barzo, bardzo, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim = c(0,1), xlab = "year", ylab = "proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.red, lwd=3)
legend("topleft", c("barzo > bardzo"), 
       text.col=my.red,
       bty="n",
       lty=1, col=my.red, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))








# Fig. 3. the change -bych/-bychmy > -bym/-byśmy
#
piotrowski = uporzadkuj(bych, bym, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.green, lwd=3)
legend("topleft", c("-bych/-bychmy > -bym/-byśmy"), 
       text.col=my.green,
       bty="n",
       lty=1, col=my.green, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))
```






# Fig. 4. The course of change of the form _bychmy_ > _byśmy_ alone.
#
piotrowski = uporzadkuj(bychmy, bysmy, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.violet, lwd=3)
legend("bottomleft", c("-bychmy > -byśmy"), 
       text.col=my.violet,
       bty="n",
       lty=1, col=my.violet, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))
```








# Fig. 5. The course of change of the superlative marker _na-_ > _naj-_.
#
piotrowski = uporzadkuj(stopniowanie_na, stopniowanie_naj, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.grey, lwd=3)
legend("topleft", c("na- > naj-"), 
       text.col=my.grey,
       bty="n",
       lty=1, col=my.grey, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))








# Fig. 6. The course of change -_ir_- > -_er_-.
#
piotrowski = uporzadkuj(ir, ier, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.violet, lwd=3)
legend("topleft", c("-ir- > -er-"), 
       text.col=my.violet,
       bty="n",
       lty=1, col=my.violet, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))









# Fig. 7. The course of change _inszy_ > _inny_ (in all inflected forms).
#
piotrowski = uporzadkuj(inszy, inny, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim = c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),3), family = quasibinomial(logit))
lines(names(piotrowski), model$fitted, type = "l", col = my.orange, lwd = 3)
legend("topleft", c("inszy > inny"), 
       text.col = my.orange,
       bty = "n",
       lty = 1, col = my.orange, lwd = 3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))









# Fig. 8. The course of change _wszytek_ > _wszystek_ (in all inflecting forms).
#
piotrowski = uporzadkuj(wszytko, wszystko, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),3), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.teal, lwd=3)
legend("topleft", c("wszytek > wszystek"), 
       text.col=my.teal,
       bty="n",
       lty=1, col=my.teal, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("bottomright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))









# Fig. 9. The competing forms _abo_ and _albo_: a model of polynomial logistic regression.
#
piotrowski = uporzadkuj(albo, abo, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),6), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.violet, lwd=3)
legend("topleft", c("abo > albo"), 
       text.col=my.violet,
       bty="n",
       lty=1, col=my.violet, lwd=3
       )
r2 = round(NagelkerkeR2(model)$R2, 3)
legend("topright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke), list(Nagelkerke = r2)))










# Fig. 10. The competing forms _abo_ and _albo_: two independent logistic models computed separately for two time spans, 1380–1610 and 1610–1850.
#
piotrowski = uporzadkuj(albo, abo, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form")

aa = piotrowski[1:22]
model_aa = glm(aa ~ as.numeric(names(aa)), family=quasibinomial(logit))
bb = piotrowski[22:46]
model_bb = glm(bb ~ as.numeric(names(bb)), family=quasibinomial(logit))
lines(names(aa), model_aa$fitted, type="l", col=my.red, lwd=3)
lines(names(bb), model_bb$fitted, type="l", col=my.teal, lwd=3)
abline(v = 1610, lty = 2)
legend("topleft", c("abo > albo (a)", "abo > albo (b)"), 
       text.col = c(my.red, my.teal),
       bty = "n",
       lty = 1, col = c(my.red, my.teal), lwd=3
       )
r2_aa = round(NagelkerkeR2(model_aa)$R2, 3)
r2_bb = round(NagelkerkeR2(model_bb)$R2, 3)
legend("topright", bty="n", legend = substitute(paste(italic(R)^2, " = ", Nagelkerke, " (a),    ", paste(italic(R)^2, " = ", Nagelkerke2), " (b)"), list(Nagelkerke = r2_aa, Nagelkerke2 = r2_bb)))












# Fig. 11. The overall picture of the language changes in Middle Polish.
#
# więtszy/większy
piotrowski = uporzadkuj(wietszy, wiekszy, zakres.chrono = 20, krok = 10)
plot(piotrowski ~ as.numeric(names(piotrowski)), ylim=c(0,1), xlab="year", ylab="proportion of the innovative form", type="n")
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.blue, lwd=3)
#
# byśmy/bychmy
piotrowski = uporzadkuj(bych, bym, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.green, lwd=3)
#
# bardzo/barzo
piotrowski = uporzadkuj(barzo, bardzo, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.red, lwd=3)
#
# na/naj
piotrowski = uporzadkuj(stopniowanie_na, stopniowanie_naj, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.grey, lwd=3)
#
# inny/inszy
piotrowski = uporzadkuj(inszy, inny, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),3), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.orange, lwd=3)
#
# wszytko/wszystko
piotrowski = uporzadkuj(wszytko, wszystko, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),3), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.teal, lwd=3)
#
# abo/albo
piotrowski = uporzadkuj(albo, abo, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)),6), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.violet, lwd=3)
#
# ir/ier
piotrowski = uporzadkuj(ir, ier, zakres.chrono = 20, krok = 10)
model = glm(piotrowski ~ as.numeric(names(piotrowski)), family=quasibinomial(logit))
lines(names(piotrowski), model$fitted, type="l", col=my.violet, lwd=3)
#
legend("topleft", c("więtszy > większy", "-bych > -bym", "barzo > bardzo", "na- > naj-", "inszy > inny", "wszytek > wszystek", "abo > albo", "-ir- > -er-"), 
       text.col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal, my.violet, my.violet),
       bty="n",
       lty=1, lwd=3, col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal, my.violet, my.violet)
       )









# Fig. 12. The influence of the subcorpus size on the goodness of fit (10-year overlap).
#


dopasowanie = function(przed, po, wielomian) {
    r2_all = c()
    for(slice in seq(10,100,5)) {
        piotrowski = uporzadkuj(przed, po, zakres.chrono = slice, krok = 10)
        # sanitizing the names, so that they don't contain any non-numeric values
        names(piotrowski) = gsub("([0-9]{4}).*", "\\1", names(piotrowski))
        model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)), wielomian), family = quasibinomial(logit))
        r2 = round(NagelkerkeR2(model)$R2, 3)
        r2_all = c(r2_all, r2)
    }
    return(r2_all)
}

#
dawniejsze = list(wietszy, bych, barzo, stopniowanie_na, inszy, wszytko)
nowsze = list(wiekszy, bym, bardzo, stopniowanie_naj, inny, wszystko)
kolor = c(my.blue, my.green, my.red, my.grey, my.orange, my.teal)
daty = seq(10,100,5)
#
plot(dopasowanie(dawniejsze[[1]], nowsze[[1]], 1), type = "n", xlim = c(0,100), ylim = c(0.3,1), ylab = substitute(paste(italic(R)^2, " (Nagelkerke)")), xlab = "subcorpus size (years)", axes = F)
axis(1)
axis(2)
box()
#
for(j in 1:length(dawniejsze) ) {
    wielomian = c(1,1,1,1,4,4)
    x = dopasowanie(dawniejsze[[j]], nowsze[[j]], wielomian[j])
    lines(x ~ daty, col = kolor[j], lwd = 3)
}
#
legend("bottomright", c("więtszy > większy", "-bych > -bym", "barzo > bardzo", "na- > naj-", "inszy > inny", "wszytek > wszystek"), 
       text.col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal),
       bty="n",
       lty=1, lwd=3, col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal)
       )
#
legend("bottomleft", "overlap = 10 years", bty = "n")







# Fig. 13. The influence of the subcorpus size on the goodness of fit (20-year overlap)."}
#

dopasowanie = function(przed, po, wielomian) {
    r2_all = c()
    for(slice in seq(20,100,5)) {
        piotrowski = uporzadkuj(przed, po, zakres.chrono = slice, krok = 20)
        # sanitizing the names, so that they don't contain any non-numeric values
        names(piotrowski) = gsub("([0-9]{4}).*", "\\1", names(piotrowski))
        model = glm(piotrowski ~ poly(as.numeric(names(piotrowski)), wielomian), family=quasibinomial(logit))
        r2 = round(NagelkerkeR2(model)$R2, 3)
        r2_all = c(r2_all, r2)
    }
    return(r2_all)
}

#
dawniejsze = list(wietszy, bych, barzo, stopniowanie_na, inszy, wszytko)
nowsze = list(wiekszy, bym, bardzo, stopniowanie_naj, inny, wszystko)
kolor = c(my.blue, my.green, my.red, my.grey, my.orange, my.teal)
daty = seq(20,100,5)
#
plot(dopasowanie(dawniejsze[[1]], nowsze[[1]], 1), type = "n", xlim = c(0,100), ylim = c(0.3,1), ylab = substitute(paste(italic(R)^2, " (Nagelkerke)")), xlab = "subcorpus size (years)", axes = F)
axis(1)
axis(2)
box()
#
for(j in 1:length(dawniejsze) ) {
    wielomian = c(1,1,1,1,4,4)
    x = dopasowanie(dawniejsze[[j]], nowsze[[j]], wielomian[j])
    lines(x ~ daty, col = kolor[j], lwd = 3)
}
#
legend("bottomright", c("więtszy > większy", "-bych > -bym", "barzo > bardzo", "na- > naj-", "inszy > inny", "wszytek > wszystek"), 
       text.col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal),
       bty="n",
       lty=1, lwd=3, col=c(my.blue, my.green, my.red, my.grey, my.orange, my.teal)
       )
#
legend("bottomleft", "overlap = 20 years", bty = "n")




