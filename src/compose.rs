pub struct BoxFn<'t,Args,Output>(pub Box<dyn Fn<Args,Output=Output>+'t>);
impl<Args,Output> FnOnce<Args> for BoxFn<'_,Args,Output> {
    type Output=Output;
    extern "rust-call" fn call_once(self, args:Args) -> Self::Output { self.0.call(args) }
}
impl<Args,Output> FnMut<Args> for BoxFn<'_,Args,Output> { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { self.0.call(args) } }
impl<Args,Output> Fn<Args> for BoxFn<'_,Args,Output> { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args) } }
impl<'t,Args,Output> BoxFn<'t,Args,Output> { pub fn new<F:Fn<Args,Output=Output>+'t>(f:F) -> Self { Self(Box::new(f)) } } // type alias hides constructor

macro_rules! unary { ([$($Op:ident $op:ident),+]) => ($(
    pub struct $Op<'t,Args,A>(pub BoxFn<'t,Args,A>);
    impl<Args,A> FnOnce<Args> for $Op<'_,Args,A> where A:std::ops::$Op {
        type Output = <A as std::ops::$Op>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args,A> FnMut<Args> for $Op<'_,Args,A> where A:std::ops::$Op { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args,A> Fn<Args> for $Op<'_,Args,A> where A:std::ops::$Op { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args).$op() } }
    impl<'t,Args:'static,A:'static> std::ops::$Op for BoxFn<'t,Args,A> where A:std::ops::$Op {
        type Output = BoxFn<'t,Args,<<Self as FnOnce<Args>>::Output as std::ops::$Op>::Output>;
        fn $op(self) -> Self::Output { BoxFn::new($Op(self)) }
    }
)+)}

macro_rules! binary { ([$($Op:ident $op:ident),+] $Uniform:ty) => (
mod binary_box {$(
    pub struct $Op<'t,Args,A,B>(pub super::BoxFn<'t,Args,A>, pub super::BoxFn<'t,Args,B>);
    impl<Args:Copy,A,B> FnOnce<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> {
        type Output = <A as std::ops::$Op<B>>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args:Copy,A,B> FnMut<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args:Copy,A,B> Fn<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B> { extern "rust-call" fn call(&self, args:Args) -> Self::Output { self.0.call(args).$op(self.1.call(args)) } }
)+}
$(
    impl<'t,Args:Copy+'t,B:'t,A:std::ops::$Op<B>+'t,G:Fn<Args,Output=B>+'t> std::ops::$Op<G> for BoxFn<'t,Args,A> {
        type Output = BoxFn<'t,Args, <binary_box::$Op<'t,Args,A,B> as FnOnce<Args>>::Output>;
        fn $op(self, b:G) -> Self::Output { BoxFn::new(binary_box::$Op(BoxFn::new(self),BoxFn::new(b))) }
    }
    impl<'t,Args:Copy,B:'t,A:std::ops::$Op<B>+'t,G:Fn<Args,Output=B>+'t> std::ops::$Op<G> for &'t BoxFn<'t,Args,A> {
        type Output = BoxFn<'t,Args, <binary_box::$Op<'t,Args,A,B> as FnOnce<Args>>::Output>;
        fn $op(self, b:G) -> Self::Output { BoxFn::new(binary_box::$Op(BoxFn::new(self),BoxFn::new(b))) }
    }
)+
mod uniform_binary_box {$(
    pub struct $Op<'t,Args,A,B>(pub A, pub super::BoxFn<'t,Args,B>);
    impl<Args,A,B> FnOnce<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy {
        type Output = <A as std::ops::$Op<B>>::Output;
        extern "rust-call" fn call_once(self, args:Args) -> Self::Output { Self::call(&self, args) }
    }
    impl<Args,A,B> FnMut<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy { extern "rust-call" fn call_mut(&mut self, args:Args) -> Self::Output { Self::call(&self, args) } }
    impl<Args,A,B> Fn<Args> for $Op<'_,Args,A,B> where A:std::ops::$Op<B>+Copy { extern "rust-call" fn call(&self, args:Args) -> Self::Output { (self.0).$op(self.1.call(args)) } }
)+}
$(
    impl<'t,Args:'t,B:'t> std::ops::$Op<BoxFn<'t,Args,B>> for $Uniform where Self:std::ops::$Op<B> {
        type Output = BoxFn<'t,Args,<Self as std::ops::$Op<B>>::Output>;
        fn $op(self, b:BoxFn<'t,Args,B>) -> Self::Output { BoxFn::new(uniform_binary_box::$Op(self,b)) }
    }
    impl<'t,Args,B> std::ops::$Op<&'t BoxFn<'t,Args,B>> for $Uniform where Self:std::ops::$Op<B> {
        type Output = BoxFn<'t,Args,<Self as std::ops::$Op<B>>::Output>;
        fn $op(self, b:&'t BoxFn<'t,Args,B>) -> Self::Output { BoxFn::new(uniform_binary_box::$Op(self,BoxFn::new(b))) }
    }
)+
)}

unary!([Neg neg]);
binary!([Add add, Sub sub, Mul mul] f32);
