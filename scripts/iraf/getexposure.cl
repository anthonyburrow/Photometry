#list of images to extract exposures
list = "stands"

s2 = "EXPTIME"
while (fscan(list, s1) !=EOF){
      printf (s1)
      hselect (s1,s2,yes)
}
