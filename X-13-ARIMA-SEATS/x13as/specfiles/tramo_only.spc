series{
  title = "Consumer Food Price Index - All India Combined"
  start = 2013.01
  span = (2013.01, 2024.08)
  data = (
    105.4 106.4 106.5 107.5 109.1 112.4 115.2 117.3 119.0 121.1 123.9 118.7
    115.6 114.8 115.7 117.4 118.8 120.5 125.4 127.5 126.4 125.8 125.3 123.4
    122.7 122.7 122.8 123.4 124.5 127.1 128.1 130.3 131.3 132.4 132.9 131.3
    131.1 129.2 129.2 131.3 133.8 137.0 138.8 138.0 136.5 136.8 135.6 133.1
    131.9 131.8 131.8 132.1 132.4 134.1 138.3 140.1 138.2 139.4 141.5 139.7
    138.1 136.1 135.5 135.8 136.5 138.0 140.1 140.5 138.9 138.2 137.8 136.0
    135.0 135.1 135.9 137.3 139.0 141.1 143.4 144.7 146.0 149.1 151.6 155.3
    153.4 149.7 147.8 153.4 151.8 153.4 156.7 157.8 161.6 165.5 166.0 160.6
    156.4 155.5 155.0 156.4 159.4 161.3 162.9 162.7 162.7 166.9 169.1 167.1
    164.9 164.6 166.9 169.4 172.1 173.8 173.8 175.1 176.7 178.6 177.0 174.1
    174.8 174.4 174.9 175.9 177.2 181.7 193.8 192.5 188.4 190.4 192.4 190.7
    189.3 189.5 189.8 191.2 192.6 198.7 204.3 203.4
  )
}

transform{
  function = auto
  print = all
}

automdl{print = all 
maxdiff = (2 1)}

outlier{types = (ls ao)
method = addone
print = all
}

estimate{
    print = all
}

regression {
  variables = (td1nolpyear)
  user = (diwali)
  aictest = (user td1nolpyear)
  start = 2013.01
  data = (
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0405 -0.0405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.6595  0.6595  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.3595  0.3595  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.6595  0.6595  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.0595  0.0595  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.6595  0.6595  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.3405 -0.3405  0.0)
   print = all
}
