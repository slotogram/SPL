#ifndef _SPL_WRAP_STRUCT_
#define _SPL_WRAP_STRUCT_

/* abstract wrapper macros */

#define head2(a,b) a
#define head3(a,b,c) a
#define head4(a,b,c,d) a
#define head5(a,b,c,d,e) a
#define head6(a,b,c,d,e,f) a
#define head7(a,b,c,d,e,f,g) a
#define head8(a,b,c,d,e,f,g,h) a
#define head9(a,b,c,d,e,f,g,h,i) a
#define head10(a,b,c,d,e,f,g,h,i,j) a
#define head11(a,b,c,d,e,f,g,h,i,j,k) a
#define head12(a,b,c,d,e,f,g,h,i,j,k,l) a
#define head13(a,b,c,d,e,f,g,h,i,j,k,l,m) a
#define head14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) a
#define head15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) a
#define head16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) a
#define head17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) a
#define head18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) a
#define head19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) a
#define head20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) a
#define head21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) a
#define head22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) a
#define head23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) a
#define head24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) a
#define head25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) a
#define head26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) a

#define tail2(a,b) b
#define tail3(a,b,c) b,c
#define tail4(a,b,c,d) b,c,d
#define tail5(a,b,c,d,e) b,c,d,e
#define tail6(a,b,c,d,e,f) b,c,d,e,f
#define tail7(a,b,c,d,e,f,g) b,c,d,e,f,g
#define tail8(a,b,c,d,e,f,g,h) b,c,d,e,f,g,h
#define tail9(a,b,c,d,e,f,g,h,i) b,c,d,e,f,g,h,i
#define tail10(a,b,c,d,e,f,g,h,i,j) b,c,d,e,f,g,h,i,j
#define tail11(a,b,c,d,e,f,g,h,i,j,k) b,c,d,e,f,g,h,i,j,k
#define tail12(a,b,c,d,e,f,g,h,i,j,k,l) b,c,d,e,f,g,h,i,j,k,l
#define tail13(a,b,c,d,e,f,g,h,i,j,k,l,m) b,c,d,e,f,g,h,i,j,k,l,m
#define tail14(a,b,c,d,e,f,g,h,i,j,k,l,m,n) b,c,d,e,f,g,h,i,j,k,l,m,n
#define tail15(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) b,c,d,e,f,g,h,i,j,k,l,m,n,o
#define tail16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p
#define tail17(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q
#define tail18(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r
#define tail19(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s
#define tail20(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t
#define tail21(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u
#define tail22(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v
#define tail23(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w
#define tail24(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x
#define tail25(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y
#define tail26(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z

#define join2(sep,args)  head2  args WRAP##sep tail2 args
#define join3(sep,args)  head3  args WRAP##sep join2 (sep,(tail3  args))
#define join4(sep,args)  head4  args WRAP##sep join3 (sep,(tail4  args))
#define join5(sep,args)  head5  args WRAP##sep join4 (sep,(tail5  args))
#define join6(sep,args)  head6  args WRAP##sep join5 (sep,(tail6  args))
#define join7(sep,args)  head7  args WRAP##sep join6 (sep,(tail7  args))
#define join8(sep,args)  head8  args WRAP##sep join7 (sep,(tail8  args))
#define join9(sep,args)  head9  args WRAP##sep join8 (sep,(tail9  args))
#define join10(sep,args) head10 args WRAP##sep join9 (sep,(tail10 args))
#define join11(sep,args) head11 args WRAP##sep join10(sep,(tail11 args))
#define join12(sep,args) head12 args WRAP##sep join11(sep,(tail12 args))
#define join13(sep,args) head13 args WRAP##sep join12(sep,(tail13 args))
#define join14(sep,args) head14 args WRAP##sep join13(sep,(tail14 args))
#define join15(sep,args) head15 args WRAP##sep join14(sep,(tail15 args))
#define join16(sep,args) head16 args WRAP##sep join15(sep,(tail16 args))
#define join17(sep,args) head17 args WRAP##sep join16(sep,(tail17 args))
#define join18(sep,args) head18 args WRAP##sep join17(sep,(tail18 args))
#define join19(sep,args) head19 args WRAP##sep join18(sep,(tail19 args))
#define join20(sep,args) head20 args WRAP##sep join19(sep,(tail20 args))
#define join21(sep,args) head21 args WRAP##sep join10(sep,(tail21 args))
#define join22(sep,args) head22 args WRAP##sep join11(sep,(tail22 args))
#define join23(sep,args) head23 args WRAP##sep join12(sep,(tail23 args))
#define join24(sep,args) head24 args WRAP##sep join13(sep,(tail24 args))
#define join25(sep,args) head25 args WRAP##sep join14(sep,(tail25 args))
#define join26(sep,args) head26 args WRAP##sep join15(sep,(tail26 args))

#define join(sep,nargs,args) join##nargs(sep,args)

/* helper macro definitions */

#define WRAP
#define WRAPCOMMA ,
#define WRAPSEMICOLON ;
#define WRAPSTRUCT(s) , (s)->
#define WRAPRETURN(type) WRAPRETURN_##type
#define WRAPRETURN_void
#define WRAPRETURN_unsigned return

/* wrapper names */

#define wrap_func_name(name) name##_str
#define wrap_struct_name(name) _##name##_args

/* declare wrapper functions */

#define declare_wrap_func(ret,name,nargs,args) \
struct wrap_struct_name(name) { join(SEMICOLON,nargs,args); }; \
EXPORT ret name args; \
EXPORT ret wrap_func_name(name) ( struct wrap_struct_name(name) * );

/* define wrapper function */

#define define_wrap_func(ret,name,nargs,args) \
EXPORT ret wrap_func_name(name) ( struct wrap_struct_name(name) *s ) { WRAPRETURN(ret) name((s)-> join(STRUCT(s),nargs,args)); }

#endif//_SPL_WRAP_STRUCT_
