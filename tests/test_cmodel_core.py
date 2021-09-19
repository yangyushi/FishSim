from scipy.spatial.transform import Rotation
from itertools import product
import fish_sim as fs
import numpy as np



def test_unrave_index():
    data_2d = np.empty((
        np.random.randint(5, 50),
        np.random.randint(5, 50)
    ))
    data_3d = np.empty((
        np.random.randint(5, 50),
        np.random.randint(5, 50),
        np.random.randint(5, 50)
    ))
    for _ in range(10):
        index_1d = np.random.randint(0, data_2d.size)
        i_cpp = fs.cmodel._unravel_index_2d(index_1d, np.array(data_2d.shape))
        i_np = np.unravel_index(index_1d, data_2d.shape)
        for i1, i2 in zip(i_cpp, i_np):
            assert i1 == i2
    print("2D unravel index okay")

    for _ in range(10):
        index_1d = np.random.randint(0, data_3d.size)
        i_cpp = fs.cmodel._unravel_index_3d(index_1d, np.array(data_3d.shape))
        i_np = np.unravel_index(index_1d, data_3d.shape)
        for i1, i2 in zip(i_cpp, i_np):
            assert i1 == i2
    print("3D unravel index okay")


def test_rotation():
    for _ in range(100):
        v1 = np.random.random(3) - 0.5
        R = fs.cmodel._get_rotation_matrix_from_001(*v1)
        v2 = R @ np.array((0, 0, 1)).reshape((3, 1))
        v2 = np.squeeze(v2)
        v2 = v2 / np.linalg.norm(v2) * np.linalg.norm(v1)
        assert np.isclose(v1, v2).all()
    print("rotation from 001 okay")

    for _ in range(100):
        v1 = np.random.random(3) - 0.5
        v2 = np.random.random(3) - 0.5
        R = fs.cmodel._get_rotation_matrix(v1, v2)
        v3 = R @ v1
        v3 = v3 / np.linalg.norm(v3) * np.linalg.norm(v2)
        assert np.isclose(v2, v3).all()
    print("rotation between two vectors okay")

    for _ in range(100):
        theta_lim = np.random.uniform(0, np.pi)
        v1 = np.random.random(3) - 0.5
        v2 = np.random.random(3) - 0.5
        angle = np.arccos((v1 @ v2) / np.linalg.norm(v1) / np.linalg.norm(v2))
        R_lim = fs.cmodel._get_rotation_matrix_limited(v1, v2, theta_lim)
        theta = min(theta_lim, angle)  # the actual angle being rotated
        R_vec = Rotation.from_matrix(R_lim).as_rotvec()
        assert np.isclose(np.linalg.norm(R_vec), theta)
        R_vec = R_vec / theta * angle
        R = Rotation.from_rotvec(R_vec).as_matrix()
        v3 = R @ v1
        v3 = v3 / np.linalg.norm(v3) * np.linalg.norm(v2)
        assert np.isclose(v2, v3).all()
    print("limited rotation between two vectors okay")


def test_product():
    for _ in range(10):
        size = np.random.randint(0, 10)
        a1 = np.random.randint(0, 100, size)
        a2 = np.random.randint(0, 100, size)
        a3 = np.random.randint(0, 100, size)

        product_py = list(product(a1, a2))
        product_cpp = fs.cmodel._product_2d(a1, a2)
        for p1, p2 in zip(product_py, product_cpp):
            for v1, v2 in zip(p1, p2):
                assert v1 == v2

        product_py = list(product(a1, a2, a3))
        product_cpp = fs.cmodel._product_3d(a1, a2, a3)
        for p1, p2 in zip(product_py, product_cpp):
            for v1, v2 in zip(p1, p2):
                assert v1 == v2

    print("Cartesian product 2D/3D okay")

if __name__ == "__main__":
    test_unrave_index()
    test_rotation()
    test_product()
