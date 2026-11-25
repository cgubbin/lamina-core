use crate::{
    ComplexScalar, Constant, DerivativeOrder, DifferentiableMaterial,
    DifferentiableMeromorphicMaterial, MeromorphicMaterial, Sampled, material::Material,
};

/// Linear dispersive material used for testing.
///
/// ε(k₀) = ε₀ + slope · k₀
/// μ(k₀) = μ₀ + slope · k₀
#[derive(Clone, Debug)]
pub struct LinearDispersion {
    pub(super) epsilon0: f64,
    pub(super) epsilon_slope: f64,
    pub(super) mu0: f64,
    pub(super) mu_slope: f64,
}

impl LinearDispersion {
    pub fn new(epsilon0: f64, epsilon_slope: f64, mu0: f64, mu_slope: f64) -> Self {
        Self {
            epsilon0,
            epsilon_slope,
            mu0,
            mu_slope,
        }
    }
}

pub fn vacuum() -> Constant<f64> {
    Constant::new(1.0, 1.0)
}

pub fn constant(epsilon: f64, mu: f64) -> Constant<f64> {
    Constant::new(epsilon, mu)
}

pub fn linear(epsilon0: f64, epsilon_slope: f64, mu0: f64, mu_slope: f64) -> LinearDispersion {
    LinearDispersion::new(epsilon0, epsilon_slope, mu0, mu_slope)
}

impl Material for LinearDispersion {
    type Real = f64;

    fn relative_permittivity<I, C>(&self, wavenumber: I) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        wavenumber.map(|x| {
            C::from_real(self.epsilon0) + C::from_real(x) * C::from_real(self.epsilon_slope)
        })
    }

    fn relative_permeability<I, C>(&self, wavenumber: I) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        wavenumber.map(|x| C::from_real(self.mu0) + C::from_real(x) * C::from_real(self.mu_slope))
    }
}

impl DifferentiableMaterial for LinearDispersion {
    fn relative_permittivity_derivative<I, C>(
        &self,
        wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        match order {
            DerivativeOrder::First => wavenumber.map(|_| C::from_real(self.epsilon_slope)),
            DerivativeOrder::Second => wavenumber.map(|_| C::zero()),
            DerivativeOrder::Third => wavenumber.map(|_| C::zero()),
        }
    }

    fn relative_permeability_derivative<I, C>(
        &self,
        wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        match order {
            DerivativeOrder::First => wavenumber.map(|_| C::from_real(self.mu_slope)),
            DerivativeOrder::Second => wavenumber.map(|_| C::zero()),
            DerivativeOrder::Third => wavenumber.map(|_| C::zero()),
        }
    }
}

impl MeromorphicMaterial for LinearDispersion {
    fn relative_permittivity_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber
            .map(|x| C::from_real(self.epsilon0) + x * C::from_real(self.epsilon_slope))
    }

    fn relative_permeability_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|x| C::from_real(self.mu0) + x * C::from_real(self.mu_slope))
    }
}

impl DifferentiableMeromorphicMaterial for LinearDispersion {
    fn relative_permittivity_complex_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        match order {
            DerivativeOrder::First => vacuum_wavenumber.map(|_| C::from_real(self.epsilon_slope)),
            DerivativeOrder::Second => vacuum_wavenumber.map(|_| C::zero()),
            DerivativeOrder::Third => vacuum_wavenumber.map(|_| C::zero()),
        }
    }

    fn relative_permeability_complex_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C>,
    {
        match order {
            DerivativeOrder::First => vacuum_wavenumber.map(|_| C::from_real(self.mu_slope)),
            DerivativeOrder::Second => vacuum_wavenumber.map(|_| C::zero()),
            DerivativeOrder::Third => vacuum_wavenumber.map(|_| C::zero()),
        }
    }
}

#[derive(Clone, Debug)]
pub struct QuadraticDispersion {
    pub(super) epsilon0: f64,
    pub(super) epsilon_slope: f64,
    pub(super) epsilon_curvature: f64,
    pub(super) mu0: f64,
    pub(super) mu_slope: f64,
    pub(super) mu_curvature: f64,
}

impl QuadraticDispersion {
    pub fn new(
        epsilon0: f64,
        epsilon_slope: f64,
        epsilon_curvature: f64,
        mu0: f64,
        mu_slope: f64,
        mu_curvature: f64,
    ) -> Self {
        Self {
            epsilon0,
            epsilon_slope,
            epsilon_curvature,
            mu0,
            mu_slope,
            mu_curvature,
        }
    }
}

pub fn quadratic(
    epsilon0: f64,
    epsilon_slope: f64,
    epsilon_curvature: f64,
    mu0: f64,
    mu_slope: f64,
    mu_curvature: f64,
) -> QuadraticDispersion {
    QuadraticDispersion::new(
        epsilon0,
        epsilon_slope,
        epsilon_curvature,
        mu0,
        mu_slope,
        mu_curvature,
    )
}

impl Material for QuadraticDispersion {
    type Real = f64;

    fn relative_permittivity<I, C>(&self, wavenumber: I) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        wavenumber.map(|x| {
            C::from_real(self.epsilon0)
                + C::from_real(x) * C::from_real(self.epsilon_slope)
                + C::from_real(x) * C::from_real(x) * C::from_real(self.epsilon_curvature)
        })
    }

    fn relative_permeability<I, C>(&self, wavenumber: I) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        wavenumber.map(|x| {
            C::from_real(self.mu0)
                + C::from_real(x) * C::from_real(self.mu_slope)
                + C::from_real(x) * C::from_real(x) * C::from_real(self.mu_curvature)
        })
    }
}

impl DifferentiableMaterial for QuadraticDispersion {
    fn relative_permittivity_derivative<I, C>(
        &self,
        wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        I: Sampled<Elem = Self::Real>,
        C: ComplexScalar<RealField = f64>,
    {
        match order {
            DerivativeOrder::First => wavenumber.map(|x| {
                C::from_real(self.epsilon_slope)
                    + C::from_real(2.0 * x) * C::from_real(self.epsilon_curvature)
            }),
            DerivativeOrder::Second => {
                wavenumber.map(|_| C::from_real(2.0) * C::from_real(self.epsilon_curvature))
            }
            DerivativeOrder::Third => wavenumber.map(|_| C::zero()),
        }
    }

    fn relative_permeability_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real>,
        I: Sampled<Elem = C::RealField>,
    {
        match order {
            DerivativeOrder::First => vacuum_wavenumber.map(|x| {
                C::from_real(self.mu_slope)
                    + C::from_real(2.0 * x) * C::from_real(self.mu_curvature)
            }),
            DerivativeOrder::Second => {
                vacuum_wavenumber.map(|_| C::from_real(2.0) * C::from_real(self.mu_curvature))
            }
            DerivativeOrder::Third => vacuum_wavenumber.map(|_| C::zero()),
        }
    }
}

impl MeromorphicMaterial for QuadraticDispersion {
    fn relative_permittivity_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|x| {
            C::from_real(self.epsilon0)
                + x * C::from_real(self.epsilon_slope)
                + x * x * C::from_real(self.epsilon_curvature)
        })
    }

    fn relative_permeability_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|x| {
            C::from_real(self.mu0)
                + x * C::from_real(self.mu_slope)
                + x * x * C::from_real(self.mu_curvature)
        })
    }
}

impl DifferentiableMeromorphicMaterial for QuadraticDispersion {
    fn relative_permittivity_complex_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        match order {
            DerivativeOrder::First => vacuum_wavenumber.map(|x| {
                C::from_real(self.epsilon_slope)
                    + C::from_real(2.0) * x * C::from_real(self.epsilon_curvature)
            }),
            DerivativeOrder::Second => {
                vacuum_wavenumber.map(|_| C::from_real(2.0) * C::from_real(self.epsilon_curvature))
            }
            DerivativeOrder::Third => vacuum_wavenumber.map(|_| C::zero()),
        }
    }

    fn relative_permeability_complex_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = Self::Real> + Copy,
        I: Sampled<Elem = C>,
    {
        match order {
            DerivativeOrder::First => vacuum_wavenumber.map(|x| {
                C::from_real(self.mu_slope)
                    + C::from_real(2.0) * x * C::from_real(self.mu_curvature)
            }),
            DerivativeOrder::Second => {
                vacuum_wavenumber.map(|_| C::from_real(2.0) * C::from_real(self.mu_curvature))
            }
            DerivativeOrder::Third => vacuum_wavenumber.map(|_| C::zero()),
        }
    }
}

/// `ε = Σₙ aₙ k₀ⁿ` and `μ = Σₙ bₙ k₀ⁿ`, `n ≤ 2`, with complex coefficients.
#[derive(Clone, Copy, Debug)]
pub struct ComplexQuadratic {
    pub epsilon: [num_complex::Complex64; 3],
    pub mu: [num_complex::Complex64; 3],
}

/// The value (`order` `None`) or a derivative of `Σₙ cₙ kⁿ` at `k`.
fn polynomial<X>(
    coefficients: [num_complex::Complex64; 3],
    k: X,
    order: Option<DerivativeOrder>,
) -> X
where
    X: ComplexScalar<RealField = f64>,
{
    let [c0, c1, c2] = coefficients.map(|c| X::from_real(c.re) + X::i() * X::from_real(c.im));
    let two = X::from_real(2.0);
    match order {
        None => c0 + c1 * k + c2 * k * k,
        Some(DerivativeOrder::First) => c1 + two * c2 * k,
        Some(DerivativeOrder::Second) => two * c2,
        Some(DerivativeOrder::Third) => X::zero(),
    }
}

macro_rules! complex_quadratic {
    ($trait:ident, $k:ty, $lift:path, { $($item:item)* } $($method:ident => $field:ident $(($order:ident))?;)+) => {
        impl $trait for ComplexQuadratic {
            $($item)*
            $(
                fn $method<I, X>(&self, k: I $(, $order: DerivativeOrder)?) -> I::Mapped<X>
                where
                    X: ComplexScalar<RealField = f64>,
                    I: Sampled<Elem = $k>,
                {
                    #[allow(unused_variables)]
                    let order: Option<DerivativeOrder> = None $(.or(Some($order)))?;
                    k.map(|k| polynomial(self.$field, $lift(k), order))
                }
            )+
        }
    };
}

complex_quadratic!(Material, f64, X::from_real, { type Real = f64; }
    relative_permittivity => epsilon;
    relative_permeability => mu;
);
complex_quadratic!(DifferentiableMaterial, f64, X::from_real, {}
    relative_permittivity_derivative => epsilon(order);
    relative_permeability_derivative => mu(order);
);
complex_quadratic!(MeromorphicMaterial, X, std::convert::identity, {}
    relative_permittivity_complex => epsilon;
    relative_permeability_complex => mu;
);
complex_quadratic!(DifferentiableMeromorphicMaterial, X, std::convert::identity, {}
    relative_permittivity_complex_derivative => epsilon(order);
    relative_permeability_complex_derivative => mu(order);
);

/// A tensor material with principal axes on its own frame: axis `i` has the
/// permittivity and permeability of the scalar material `self.0[i]`.
#[derive(Clone, Debug)]
pub struct Axes<M>(pub [M; 3]);

macro_rules! axes {
    ($trait:ident, $k:ty { $($item:item)* } $($method:ident => $scalar:ident $(($order:ident))?;)+) => {
        impl<M: DifferentiableMeromorphicMaterial<Real = f64>> crate::$trait for Axes<M> {
            $($item)*
            $(
                fn $method<I, X>(&self, k: I $(, $order: DerivativeOrder)?) -> I::Mapped<crate::material::Tensor3<X>>
                where
                    X: ComplexScalar<RealField = f64>,
                    I: Sampled<Elem = $k>,
                {
                    k.map(|k| {
                        let [x, y, z] = self.0.each_ref().map(|axis| {
                            axis.$scalar(crate::Scalar(k) $(, $order)?)
                        });
                        let zero = X::zero();
                        [[x, zero, zero], [zero, y, zero], [zero, zero, z]]
                    })
                }
            )+
        }
    };
}

axes!(TensorMaterial, f64 { type Real = f64; const PRINCIPAL_AXES: bool = true; }
    relative_permittivity_tensor => relative_permittivity;
    relative_permeability_tensor => relative_permeability;
);
axes!(DifferentiableTensorMaterial, f64 {}
    relative_permittivity_tensor_derivative => relative_permittivity_derivative(order);
    relative_permeability_tensor_derivative => relative_permeability_derivative(order);
);
axes!(MeromorphicTensorMaterial, X {}
    relative_permittivity_tensor_complex => relative_permittivity_complex;
    relative_permeability_tensor_complex => relative_permeability_complex;
);
axes!(DifferentiableMeromorphicTensorMaterial, X {}
    relative_permittivity_tensor_complex_derivative => relative_permittivity_complex_derivative(order);
    relative_permeability_tensor_complex_derivative => relative_permeability_complex_derivative(order);
);

/// A region of the principal-axis material `axes` oriented by `rotation`.
pub fn principal_region<M>(
    axes: [M; 3],
    rotation: crate::Rotation<f64>,
) -> crate::AnalyticalMaterialRegion<num_complex::Complex64>
where
    M: DifferentiableMeromorphicMaterial<Real = f64> + Send + Sync + 'static,
{
    crate::Constitutive::oriented(
        crate::AnalyticalTensorMaterialHandle::new(Axes(axes)),
        rotation,
    )
}
