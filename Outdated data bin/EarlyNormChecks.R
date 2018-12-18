data=read.csv(file="E:/Github_Projects/Tails/Outdated data bin/revamped_dataset.csv")
Mars=data[3:140,]
Rod=data[141:203,]

model=lm(log10( Mars$Tlength_mm)~log10(Mars$Blength_mm))
plot(model)#looks good except obviously teh wombats

model2=lm(log10( Mars$Tlength_mm)~log10(Mars$Blength_mm)*Mars$Locomotor_mode)



qqnorm(log10(Mars$Tlength_mm));qqline(log10(Mars$Tlength_mm))
hist(log10(Mars$Tlength_mm))
#taking out bandicoots
hist(log10(Mars$Tlength_mm[-(grep("Macrotis_lagotis", Mars$Species):grep("Echymipera_rufescens", Mars$Species))]))


plot(log10( Mars$Tlength_mm)~log10(Mars$Blength_mm))#; text(log10( Mars$Tlength_mm)~log10(Mars$Blength_mm), labels=Mars$Species)



Model=gls(log10 (Tlength_mm)~log10(Blength_mm), Mars)

qqnorm(log10(Mars$Blength_mm));qqline(log10(Mars$Blength_mm))
hist(log(Mars$Blength_mm))

qqnorm(log10(Rod$Tlength_mm));qqline(log10(Rod$Tlength_mm))
hist(log(Rod$Tlength_mm))

qqnorm(log10(Rod$Blength_mm));qqline(log10(Rod$Blength_mm))
hist(log(Rod$Blength_mm))
