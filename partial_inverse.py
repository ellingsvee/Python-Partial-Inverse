import scipy.sparse as sparse
import partial_inverse_module as pim


def pinv(
    M: sparse.csc_matrix, handle_singular: bool = True, raise_error: bool = True
) -> sparse.csc_matrix:
    """
    Compute the sparse partial inverse of a matrix.

    Parameters
    ----------
    M : sparse.csc_matrix
        The matrix to invert.
    handle_singular : bool, optional
        Whether to handle singular matrices, by default True.
    raise_error : bool, optional
        Whether to raise an error if the partial inverse cannot be computed, by default True.

    Returns
    -------
    sparse.csc_matrix
        The sparse partial inverse of Q.
    """

    try:
        Minv = pim.pinv(M)
    except:
        Id = sparse.eye(M.shape[0], dtype=M.dtype, format="csc")
        if handle_singular:
            M += 1e-6 * Id
            try:
                Minv = pim.pinv(M)
            except:
                if raise_error:
                    raise ValueError(
                        "handle_singular: Unable to compute the partial inverse!"
                    )
                else:
                    return Id
        else:
            if raise_error:
                raise ValueError("Matix is singular!")
            else:
                return Id
    return Minv
