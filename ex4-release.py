import decimal
a=None
w=range
g=True
O=abs
i=print
z=0.0001
def S(f,x,h=a):
 if(h==a):
  h=z
 return(f(x+h)-f(x-h))/(2*h)
def G(f,x,h=a):
 if(h==a):
  h=z
 return(f(x+h)-2*f(x)+f(x-h))/(h*h)
def N(f,x,level=1,h=a):
 if(level==0):
  return f(x)
 elif(level==2):
  return G(f,x,h)
 elif(level==1):
  return S(f,x,h)
 return N(lambda x:S(f,x),x,level-1,h)
W={0:1,1:1}
def F(n,isCache=a):
 if(n==0):
  return 1
 if(n>500 and isCache==a):
  for i in w(249,n//2):
   F(i*2,g)
 l=W.get(n,a)
 if(l==a):
  W[n]=n*F(n-1,isCache)
 return W[n]
M={0:1,1:1,2:2}
def y(n,acc=1,isCache=a):
 if(n==0 or n==1):
  return 1
 elif(n==2):
  return 2
 if(acc!=1):
  while(n!=0 and n!=1):
   acc=acc*n
   n=n-2
  return acc
 if(n>250 and isCache==a):
  for i in w(124,n//2):
   y(i*2+n%2,isCache=g)
 l=M.get(n,a)
 if(l==a):
  M[n]=n*y(n-2,isCache=isCache)
 return M[n]
def U(a,b,acc=1):
 H=acc
 if(a==b):
  return a
 if(a<b):
  for i in w(a+1,b+1):
   H=H*i
 elif(a>b):
  for i in w(b+1,a+1):
   H=H/i
 return H
m={0:{0:1},1:{0:0}}
def r(n,x,isCache=a):
 global m
 if(n>500 and isCache==a):
  for i in w(249,n//2):
   r(i*2,x,g)
 if(n==0):
  return 1
 elif(n==1):
  return x
 l=m.get(n,a)
 if(l==a):
  m[n]={}
  m[n][x]=(2*n-1)/n*x*r(n-1,x,isCache)-(n-1)/n*r(n-2,x,isCache)
 else:
  l=m[n].get(x,a)
  if(l==a):
   m[n][x]=(2*n-1)/n*x*r(n-1,x,isCache)-(n-1)/n*r(n-2,x,isCache)
 return m[n][x]
j={}
def R(l,m,x,isCache=a):
 global j
 if(O(m)>l):
  return 0
 if(m==0):
  return r(l,x)
 elif(m<0):
  return(-1)**m*U(l-m,l+m,R(l,(-1)*m,x,isCache=isCache))
 if(l==m):
  return(-1)**l*y(2*l-1)*((1-x**2)**(l/2))
 if(m==1):
  if(l==1):
   return(-1)*(1-x*x)**(1/2)
  if(l==2):
   return(-1)*3*x*(1-x*x)**(1/2)
 if(l>25 and isCache==a):
  for i in w(24,l):
   R(i,m,x,g)
 E="dL_M"
 l=j.get(l,a)
 if(l==a):
  j[l]={}
  j[l][m]={}
  if(E=="dL_M"):
   j[l][m][x]=(2*l-1)/(l-m)*x*R(l-1,m,x,isCache)-(l+m-1)/(l-m)*R(l-2,m,x,isCache)
  elif(E=="dLiM"):
   j[l][m][x]=((1-x**2)**(1/2))*(-1)/(2*m)*(R(l-1,m+1,x)+(l+m-1)*(l+m)*R(l-1,m-1,x))
 else:
  l=j[l].get(m,a)
  if(l==a):
   j[l][m]={}
   if(E=="dL_M"):
    j[l][m][x]=(2*l-1)/(l-m)*x*R(l-1,m,x,isCache)-(l+m-1)/(l-m)*R(l-2,m,x,isCache)
   elif(E=="dLiM"):
    j[l][m][x]=((1-x**2)**(1/2))*(-1)/(2*m)*(R(l-1,m+1,x)+(l+m-1)*(l+m)*R(l-1,m-1,x))
  else:
   l=j[l][m].get(x,a)
   if(l==a):
    if(E=="dL_M"):
     j[l][m][x]=(2*l-1)/(l-m)*x*R(l-1,m,x,isCache)-(l+m-1)/(l-m)*R(l-2,m,x,isCache)
    elif(E=="dLiM"):
     j[l][m][x]=((1-x**2)**(1/2))*(-1)/(2*m)*(R(l-1,m+1,x)+(l+m-1)*(l+m)*R(l-1,m-1,x))
 return j[l][m][x]
pi=3.141592653589793
def t(v):
 H=1
 v =O(v)
 while(v>2*pi):
  v+=(-1)*2*pi
 if(O(v-pi/2)<1e-11 or O(v-3*pi/2)<1e-11):
  return 0
 elif(O(v-pi)<1e-11):
  return(-1)
 elif(O(v)<1e-11 or O(v-2*pi)<1e-11):
  return 1
 for i in w(1,64):
  H+=(-1)**i*v**(i*2)/F(i*2)
 return H
def b(v):
 H=0
 if(v<0):
  while(v<(-1)*2*pi):
   v+=2*pi 
 else:
  while(v>2*pi):
   v+=(-1)*2*pi
 for i in w(1,64):
  H+=(-1)**(i+1)*v**(i*2-1)/F(i*2-1)
 return H
def X(z):
 if(z.imag!=0 and z.real==0):
  return t(z.imag)+1j*b(z.imag)
 else:
  H=1
  for i in w(1,32):
   H+=z**i/F(i)
  return H
def d(l):
 H=1
 B=2*l-3
 L=2*l-1
 while(B!=-1):
  H=H*B/L 
  if(L%2==0):
   B=B-2
  L=L-1
 return H
def c(l,m,v,phi):
 if(m==l-1):
  x=t(v)
  k=(2*l+1)/(4*pi)*d(l)*x*x*(2*l-1)**2*(1-x**2)**(l-1)
  return(-1)**(l-1)*(-1 if x<0 else 1)*(k)**(1/2)*X(1j*m*phi)
 V=R(l,m,t(v))
 if(V<0):
  o=((-1)*V)**(1/2)*1j
 else:
  o=V**(1/2)
 if(O(V)>1e+150):
  k=U(l+m,l-m,(2*l+1)*V)/(4*pi)
  if(k<0):
   J=((-1)*k)**(1/2)*1j
  else:
   J=k**(1/2)
  return J*o*X(1j*m*phi)
 else:
  k=U(l+m,l-m,(2*l+1)*V**2)/(4*pi)
  return(k)**(1/2)*(1 if V>0 else(-1))*X(1j*m*phi)
l=1000
i(c(l,l-1,pi/1000,pi/6))
i(c(l,l-1,3*pi/10,pi/6))
i(c(l,l-1,501*pi/1000,pi/6))
# Created by pyminifier (https://github.com/liftoff/pyminifier)
