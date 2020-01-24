pub struct BoxFn<'a,Args,Output>(pub Box<dyn Fn<Args,Output=Output>+'a>);
impl<Args,Output> FnOnce<Args> for BoxFn<'_,Args,Output> {
    type Output=Output;
    extern "rust-call" fn call_once(self, args:Args) -> Self::Output { self.0.call(args) }
}
impl<Args,Output> FnMut<Args> for BoxFn<'_,Args,Output> { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { self.0.call(args) } }
impl<Args,Output> Fn<Args> for BoxFn<'_,Args,Output> { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args) } }
impl<'a,Args,Output> BoxFn<'a,Args,Output> { pub fn new<F:Fn<Args,Output=Output>+'a>(f:F) -> Self { Self(Box::new(f)) } } // type alias hides constructor

macro_rules! unary { ([$($Op:ident $op:ident),+]) => ($(
    pub struct $Op<'a,Args,A>(pub super::BoxFn<'a,Args,A>);
    impl<Args,A> FnOnce<Args> for $Op<'_,Args,A> where A:std::ops::$Op {
        type Output = <A as std::ops::$Op>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args,A> FnMut<Args> for $Op<'_,Args,A> where A:std::ops::$Op { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args,A> Fn<Args> for $Op<'_,Args,A> where A:std::ops::$Op { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args).$op() } }
    impl<'a,Args:'static,A:'static> std::ops::$Op for BoxFn<'a,Args,A> where A:std::ops::$Op {
        type Output = BoxFn<'a,Args,<<Self as FnOnce<Args>>::Output as std::ops::$Op>::Output>;
        fn $op(self) -> Self::Output { BoxFn::new($Op(self)) }
    }
)+)}

macro_rules! binary { ([$($Op:ident $op:ident),+] $Uniform:ty) => (
mod binary_box {$(
    pub struct $Op<'a,Args,A,B>(pub super::BoxFn<'a,Args,A>, pub super::BoxFn<'a,Args,B>);
    impl<Args:Copy,A,B> FnOnce<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> {
        type Output = <A as std::ops::$Op<B>>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args:Copy,A,B> FnMut<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args:Copy,A,B> Fn<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args).$op(self.1.call(args)) } }
)+}
$(
    impl<'a,Args:Copy+'a,B:'a,A:std::ops::$Op<B>+'a,G:Fn<Args,Output=B>+'a> std::ops::$Op<G> for BoxFn<'a,Args,A> {
        type Output = BoxFn<'a,Args, <binary_box::$Op<'a,Args,A,B> as FnOnce<Args>>::Output>;
        fn $op(self, b:G) -> Self::Output { BoxFn::new(binary_box::$Op(BoxFn::new(self),BoxFn::new(b))) }
    }
    impl<'a,Args:Copy,B:'a,A:std::ops::$Op<B>+'a,G:Fn<Args,Output=B>+'a> std::ops::$Op<G> for &'a BoxFn<'a,Args,A> {
        type Output = BoxFn<'a,Args, <binary_box::$Op<'a,Args,A,B> as FnOnce<Args>>::Output>;
        fn $op(self, b:G) -> Self::Output { BoxFn::new(binary_box::$Op(BoxFn::new(self),BoxFn::new(b))) }
    }
)+
mod uniform_binary_box {$(
    pub struct $Op<'a,Args,A,B>(pub A, pub super::BoxFn<'a,Args,B>);
    impl<Args,A,B> FnOnce<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy {
        type Output = <A as std::ops::$Op<B>>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args,A,B> FnMut<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args,A,B> Fn<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy { extern "rust-call" fn call(&self, args:Args) -> Self::Output { (self.0).$op(self.1.call(args)) } }
)+}
$(
    impl<'a,Args:'a,B:'a> std::ops::$Op<BoxFn<'a,Args,B>> for $Uniform where Self:std::ops::$Op<B> {
        type Output = BoxFn<'a,Args,<Self as std::ops::$Op<B>>::Output>;
        fn $op(self, b:BoxFn<'a,Args,B>) -> Self::Output { BoxFn::new(uniform_binary_box::$Op(self,b)) }
    }
    impl<'a,Args,B> std::ops::$Op<&'a BoxFn<'a,Args,B>> for $Uniform where Self:std::ops::$Op<B> {
        type Output = BoxFn<'a,Args,<Self as std::ops::$Op<B>>::Output>;
        fn $op(self, b:&'a BoxFn<'a,Args,B>) -> Self::Output { BoxFn::new(uniform_binary_box::$Op(self,BoxFn::new(b))) }
    }
)+
)}

unary!([Neg neg]);
binary!([Add add, Sub sub, Mul mul] f32);
