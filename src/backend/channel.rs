//! The field subspace a backend solves.

use thiserror::Error;

use crate::{Polarisation, stack::StackSide};

/// The field subspace a backend solves: one scalar polarisation, or the
/// coupled TE/TM field of a stack whose regions mix them.
///
/// TE and TM are invariant subspaces of a *separable* stack, one whose
/// stack-frame tensors have no `xy`, `yx`, `yz` or `zy` component. The
/// coupled channel carries the four tangential components
/// `(E_x, Z₀H_y, E_y, −Z₀H_x)` of any stack.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Channel {
    /// The TE field `(E_y, H_x, H_z)`.
    TransverseElectric,
    /// The TM field `(H_y, E_x, E_z)`.
    TransverseMagnetic,
    /// The coupled field of all six components.
    Coupled,
}

impl Channel {
    /// The polarisation of a scalar channel, or `None` for the coupled one.
    pub const fn polarisation(self) -> Option<Polarisation> {
        match self {
            Self::TransverseElectric => Some(Polarisation::TransverseElectric),
            Self::TransverseMagnetic => Some(Polarisation::TransverseMagnetic),
            Self::Coupled => None,
        }
    }
}

impl From<Polarisation> for Channel {
    fn from(polarisation: Polarisation) -> Self {
        match polarisation {
            Polarisation::TransverseElectric => Self::TransverseElectric,
            Polarisation::TransverseMagnetic => Self::TransverseMagnetic,
        }
    }
}

/// A channel the backend or the stack cannot solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum ChannelError {
    /// The backend does not solve this channel.
    #[error("the {channel:?} channel is not supported by this backend")]
    Unsupported { channel: Channel },

    /// An exterior supplies a scalar branch for the coupled channel, or a pair
    /// for a scalar one.
    #[error("the {side:?} exterior branch does not match the {channel:?} channel")]
    BranchKind { side: StackSide, channel: Channel },
}
