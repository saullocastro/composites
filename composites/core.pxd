cimport numpy as np

ctypedef np.int64_t cINT
ctypedef np.double_t cDOUBLE

cdef class LaminationParameters(object):
    cdef public double xiA1, xiA2, xiA3, xiA4
    cdef public double xiB1, xiB2, xiB3, xiB4
    cdef public double xiD1, xiD2, xiD3, xiD4
    cdef public double xiE1, xiE2, xiE3, xiE4

cdef class Laminate(object):
    cdef public double A11, A12, A16, A22, A26, A66
    cdef public double B11, B12, B16, B22, B26, B66
    cdef public double D11, D12, D16, D22, D26, D66
    cdef public double E44, E45, E55
    cdef public double e1, e2, g12, nu12, nu21
    cdef public double scf_k13, scf_k23, h, offset, rho, intrho, intrhoz, intrhoz2
    cdef public list plies
    cdef public list stack
    cdef cDOUBLE[:, :] get_A(Laminate)
    cdef cDOUBLE[:, :] get_B(Laminate)
    cdef cDOUBLE[:, :] get_D(Laminate)
    cdef cDOUBLE[:, :] get_E(Laminate)
    cdef cDOUBLE[:, :] get_ABD(Laminate)
    cdef cDOUBLE[:, :] get_ABDE(Laminate)
    cpdef void rebuild(Laminate)
    cpdef void calc_scf(Laminate)
    cpdef void calc_equivalent_properties(Laminate)
    cpdef void calc_constitutive_matrix(Laminate)
    cpdef void force_orthotropic(Laminate)
    cpdef void force_symmetric(Laminate)
    cpdef LaminationParameters calc_lamination_parameters(Laminate)

