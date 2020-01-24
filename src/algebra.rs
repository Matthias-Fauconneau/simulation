use std::ops::Mul; use framework::core::Zero;

pub type Len = usize;
#[derive(Clone)] pub struct Array<T=f32, const N:Len>([T; N]);
impl<T,const N:Len> std::ops::Deref for Array<T,N> { type Target = [T; N]; fn deref(&self) -> &Self::Target { &self.0 } }
impl<T,const N:Len> std::ops::DerefMut for Array<T,N> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 } }

pub type Idx = usize;
impl<T:Copy,const N:Len> FnOnce<(Idx,)> for Array<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> FnMut<(Idx,)> for Array<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy,const N:Len> Fn<(Idx,)> for Array<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self[args.0] } }

pub fn collect<T, F:Fn(Idx)->T, const N:Len>(f : F) -> Array<T, N> { Array(framework::core::array::collect(f)) }
impl<T:Zero, const N : Idx> Zero for Array<T,N> { fn zero() -> Self { collect(|_|Zero::zero()) } }

impl<'a,T:std::ops::AddAssign<&'a T>,const N:Len> std::ops::AddAssign<&'a Array<T,N>> for Array<T,N> {
    fn add_assign(&mut self, other: &'a Self) { for (a, b) in self.iter_mut().zip(other.iter()) { *a += b; } }
}

/// Vector
pub trait Vector<T=f32,const N:Len> = Fn(Idx)->T;
pub type BoxVector<'a,T,const N:Len> = crate::compose::BoxFn<'a,(Idx,),T>;
// impl scalar · vector for Array
//pub fn mul<'a, T:Copy, const N:Len>(a:f32, b:&'a Array<T,N>) -> BoxVector<'a,T,N> where f32:Mul<T,Output=T> { BoxVector::new(move |i| a*b[i]) }

/// Operator
//#[allow(non_camel_case_types)] pub type operator<T=f32, const N:Len> = fn(&dyn Vector<T,N>)->BoxVector<T,N>;
pub trait Operator<T=f32, const N:Len> = Fn(&dyn Vector<T,N>)->BoxVector<T,N>;
pub trait OperatorOnce<T=f32, const N:Len> = FnOnce(&dyn Vector<T,N>)->BoxVector<T,N>;
pub type BoxOperator<'a, T=f32, const N:Len> = crate::compose::BoxFn<'a, (&'a dyn Vector<T,N>,), BoxVector<'a,T,N>>;

// vector · operator
impl<'a, T, const N:Len> Mul<&'a dyn Operator<T,N>> for &'a Array<f32,N> where f32:Mul<T,Output=T> {
    type Output = BoxOperator<'a,T,N>;
    fn mul(self, B:&'a dyn Operator<T,N>) -> Self::Output { BoxOperator::new(move |v| { let b=B(v); BoxVector::new(move |i| self[i] * b(i)) }) }
}
/*impl<'a, T, B:Operator<T,N>, const N:Len> Mul<B> for &'a Array<f32,N> where f32:Mul<T,Output=T> {
    type Output = BoxOperator<'a,T,N>;
    //fn mul(self, B:&'a dyn Operator<T,N>) -> Self::Output { BoxOperator::new(move |v| { let b=B(v); BoxVector::new(move |i| self[i] * b(i)) }) }
    fn mul(self, B:B) -> Self::Output { BoxOperator::new(move |v| { let b=B(v); BoxVector::new(move |i| self[i] * b(i)) }) }
}*/


pub trait Dot<T> { fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T; }
impl<T:std::iter::Sum<<f32 as Mul<T>>::Output>> Dot<T> for std::ops::Range<Idx> where f32:Mul<T> {
    fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T { self.map(|k| a(k) * b(k)).sum() }
}
impl<T:std::iter::Sum<<f32 as Mul<T>>::Output>> Dot<T> for std::ops::RangeTo<Idx> where f32:Mul<T> {
    fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized>(self, a : &A, b : &B) -> T { (0..self.end).dot(a,b) }
}
pub trait DotN<T> { fn dot<A:Fn(Idx)->f32, B:Fn(Idx)->T+?Sized,const N:Len>(self, a : &A, b : &B) -> T; }
impl<T> DotN<T> for std::ops::RangeFrom<Idx> where std::ops::Range<Idx>:Dot<T> {
    fn dot<A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(self, a : &A, b : &B) -> T { (self.start..N).dot(a,b) }
}
impl<T> DotN<T> for std::ops::RangeFull where std::ops::RangeFrom<Idx>:DotN<T> {
    fn dot<A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(self, a : &A, b : &B) -> T { (0..).dot(a,b) }
}
pub fn dot<T,A:Vector<f32,N>, B:Vector<T,N>+?Sized, const N:Len>(a : &A, b : &B) -> T where std::ops::RangeFull:DotN<T> { (..).dot(a,b) }

pub trait Matrix<const N:Len> = Fn(Idx,Idx)->f32;
pub fn matrix_mul<'a, A:Matrix<N>+'a,T:std::iter::Sum,B:Vector<T,N>+?Sized,const N:Len>(a:A, b:&'a B) -> impl Vector<T,N> + 'a where f32:Mul<T,Output=T> {
    move |i| dot(&|j| a(i, j), b)
}

pub struct SparseVector<T=f32,const N : Idx> {
    indices : Vec<Idx>,
    values : Vec<T>
}
// derive[Default] req T:Default for no reason
impl<T,const N : Idx> Default for SparseVector<T,N> { fn default() -> Self { Self{indices:Default::default(),values:Default::default()} } }

impl<T,const N : Idx> SparseVector<T,N> {
    fn get(&self, i : Idx) -> Option<T> where T:Copy { Some( self.values[ self.indices.binary_search(&i).ok()? ] ) }
    fn insert(&mut self, i : Idx, v : T) -> Option<T> {
        assert!(i<N);
        match self.indices.binary_search(&i) {
            Ok(index) => Some(std::mem::replace(&mut self.values[index], v)),
            Err(index) => { self.indices.insert(index, i); self.values.insert(index, v); None }
        }
    }
}

impl<T:Copy+Zero,const N:Len> FnOnce<(Idx,)> for SparseVector<T,N> { type Output=T;  extern "rust-call" fn call_once(self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> FnMut<(Idx,)> for SparseVector<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> Fn<(Idx,)> for SparseVector<T,N> { extern "rust-call" fn call(&self, args: (Idx,)) -> Self::Output { self.get(args.0).unwrap_or(T::zero()) } }

pub struct CSC<T=f32,const N : Idx> { columns : [SparseVector<T,N>; N as usize] }
impl<T,const N : Idx> Zero for CSC<T,N> { fn zero() -> Self { Self{columns: framework::core::array::collect(|_|Default::default())} } }

impl<T:Copy+Zero,const N : Idx> CSC<T,N> {
    fn column(&self, j : Idx) -> &SparseVector<T,N> { &self.columns[j as usize] }
    fn column_mut(&mut self, j : Idx) -> &mut SparseVector<T,N> { &mut self.columns[j as usize] }
    fn line(&self, i : Idx) -> impl Vector<T,N> + '_ { move |j| self(i, j) }
    fn get(&self, i : Idx, j : Idx) -> Option<T> { self.column(j).get(i) }

    fn insert(&mut self, i:Idx, j:Idx, v:T) -> Option<T> { self.column_mut(j).insert(i, v) }
    fn insert_once(&mut self, i:Idx, j:Idx, v:T) where T:std::fmt::Debug { self.insert(i, j, v).unwrap_none() }
}

impl<T:Copy+Zero,const N:Len> FnOnce<(Idx,Idx)> for CSC<T,N> { type Output=T; extern "rust-call" fn call_once(self, args: (Idx,Idx)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> FnMut<(Idx,Idx)> for CSC<T,N> { extern "rust-call" fn call_mut(&mut self, args: (Idx,Idx)) -> Self::Output { self.call(args) } }
impl<T:Copy+Zero,const N:Len> Fn<(Idx,Idx)> for CSC<T,N> { extern "rust-call" fn call(&self, args: (Idx,Idx)) -> Self::Output { self.get(args.0, args.1).unwrap_or(T::zero()) } }

pub type LU<const N:Len> = CSC<f32,N>;
impl<const N:Len> LU<N> {
    pub fn new<M:Matrix<N>>(A : M) -> Self {
        let mut LU = CSC::zero();
        for i in 0..N {
            for j in i..N { let v = A(i, j) - (..i).dot(&LU.line(i), LU.column(j)); LU.insert_once(i, j, v);  }
            for j in i+1..N { let v = (1. / LU(i, i)) * (A(j, i) - (..i).dot(&LU.line(i), LU.column(j))); LU.insert_once(j, i, v); }
        }
        LU // L+U-I
    }
    pub fn solve<T:Zero+Copy+std::iter::Sum+std::ops::Sub<Output=T>,B:Vector<T,N_>, const N_:Idx>(&self, b:B) -> Array<T,N_> where f32:Mul<T,Output=T> {
        assert_eq!(N, N_);
        let mut y : Array<T,N> = Zero::zero(); // Ly = b (uninit ok)
        for i in 0..N { y[i] = b(i) - (..i).dot(&self.line(i), &y); }
        let mut x : Array<T,N_> = Zero::zero(); // Ux = y (uninit ok)
        for i in N-1..=0 { x[i] = (1. / self(i, i)) * (y[i] - (i+1..N).dot(&self.line(i), &x)); }
        x
    }
}
