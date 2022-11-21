//
// Created by Kenneth Balslev on 16/11/2022.
//

#include <linalg/Matrix.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>

TEST_CASE("Matrix Tests", "[linalg]") {

    using namespace numerix::linalg;

    Matrix<int> m1(4, 4);

    m1(0,0) = 1;
    m1(0,1) = 2;
    m1(0,2) = 3;
    m1(0,3) = 4;

    m1(1,0) = 5;
    m1(1,1) = 6;
    m1(1,2) = 7;
    m1(1,3) = 8;

    m1(2,0) = 9;
    m1(2,1) = 10;
    m1(2,2) = 11;
    m1(2,3) = 12;

    m1(3,0) = 13;
    m1(3,1) = 14;
    m1(3,2) = 15;
    m1(3,3) = 16;

    SECTION("Read individual Matrix items") {

        auto m2 = m1;

        REQUIRE(m2(0,0) == 1);
        REQUIRE(m2(0,1) == 2);
        REQUIRE(m2(0,2) == 3);
        REQUIRE(m2(0,3) == 4);
        REQUIRE(m2(1,0) == 5);
        REQUIRE(m2(1,1) == 6);
        REQUIRE(m2(1,2) == 7);
        REQUIRE(m2(1,3) == 8);
        REQUIRE(m2(2,0) == 9);
        REQUIRE(m2(2,1) == 10);
        REQUIRE(m2(2,2) == 11);
        REQUIRE(m2(2,3) == 12);
        REQUIRE(m2(3,0) == 13);
        REQUIRE(m2(3,1) == 14);
        REQUIRE(m2(3,2) == 15);
        REQUIRE(m2(3,3) == 16);
    }

    SECTION("Write individual Matrix items") {

        auto m2 = m1;

        m2(0,0) = 101;
        REQUIRE(m2(0,0) == 101);

        m2(0,1) = 102;
        REQUIRE(m2(0,1) == 102);

        m2(0,2) = 103;
        REQUIRE(m2(0,2) == 103);

        m2(0,3) = 104;
        REQUIRE(m2(0,3) == 104);

        m2(1,0) = 105;
        REQUIRE(m2(1,0) == 105);

        m2(1,1) = 106;
        REQUIRE(m2(1,1) == 106);

        m2(1,2) = 107;
        REQUIRE(m2(1,2) == 107);

        m2(1,3) = 108;
        REQUIRE(m2(1,3) == 108);

        m2(2,0) = 109;
        REQUIRE(m2(2,0) == 109);

        m2(2,1) = 1010;
        REQUIRE(m2(2,1) == 1010);

        m2(2,2) = 1011;
        REQUIRE(m2(2,2) == 1011);

        m2(2,3) = 1012;
        REQUIRE(m2(2,3) == 1012);

        m2(3,0) = 1013;
        REQUIRE(m2(3,0) == 1013);

        m2(3,1) = 1014;
        REQUIRE(m2(3,1) == 1014);

        m2(3,2) = 1015;
        REQUIRE(m2(3,2) == 1015);

        m2(3,3) = 1016;
        REQUIRE(m2(3,3) == 1016);
    }

    SECTION("Read access to Matrix via iterator") {
        auto vec = std::vector<int> {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
        REQUIRE(std::equal(m1.begin(), m1.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m1.begin(), m1.end(), vec.begin()));
    }

    SECTION("Write access to Matrix via iterators") {
        auto m2 = m1;
        auto vec = std::vector<int> {101,102,103,104,105,106,107,108,109,1010,1011,1012,1013,1014,1015,1016};
        std::copy(vec.begin(), vec.end(), m2.begin());
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));
    }

    SECTION("Read individual MatrixProxy (column) items") {
        auto m2 = m1;
        auto m3 = m2({0,4,1}, {0,1,1});

        REQUIRE(m3(0,0) == 1);
        REQUIRE(m3(1,0) == 5);
        REQUIRE(m3(2,0) == 9);
        REQUIRE(m3(3,0) == 13);
    }

    SECTION("Write individual MatrixProxy (column) items") {
        auto m2 = m1;
        auto m3 = m2({0,4,1}, {0,1,1});

        m3(0,0) = 101;
        REQUIRE(m3(0,0) == 101);

        m3(1,0) = 105;
        REQUIRE(m3(1,0) == 105);

        m3(2,0) = 109;
        REQUIRE(m3(2,0) == 109);

        m3(3,0) = 1013;
        REQUIRE(m3(3,0) == 1013);

        auto vec = std::vector<int> {101,2,3,4,105,6,7,8,109,10,11,12,1013,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));
    }

    SECTION("Read access to MatrixProxy (column) via iterators") {
        auto m2 = m1;
        auto m3 = m2({0,4,1}, {0,1,1});
        auto vec = std::vector<int> {1,5,9,13};
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m3.begin(), m3.end(), vec.begin()));
    }

    SECTION("Write access to MatrixProxy (column) via iterators") {
        auto m2 = m1;
        auto m3 = m2({0,4,1}, {0,1,1});
        auto vec = std::vector<int> {101,105,109,1013};
        std::copy(vec.begin(), vec.end(), m3.begin());
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        auto vec2 = std::vector<int> {101,2,3,4,105,6,7,8,109,10,11,12,1013,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));
    }

    SECTION("Read individual MatrixProxy (row) items") {
        auto m2 = m1;
        auto m3 = m2({0,1,1}, {0,4,1});

        REQUIRE(m3(0,0) == 1);
        REQUIRE(m3(0,1) == 2);
        REQUIRE(m3(0,2) == 3);
        REQUIRE(m3(0,3) == 4);
    }

    SECTION("Write individual MatrixProxy (row) items") {
        auto m2 = m1;
        auto m3 = m2({0,1,1}, {0,4,1});

        m3(0,0) = 101;
        REQUIRE(m3(0,0) == 101);

        m3(0,1) = 102;
        REQUIRE(m3(0,1) == 102);

        m3(0,2) = 103;
        REQUIRE(m3(0,2) == 103);

        m3(0,3) = 104;
        REQUIRE(m3(0,3) == 104);

        auto vec = std::vector<int> {101,102,103,104,5,6,7,8,9,10,11,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));
    }

    SECTION("Read access to MatrixProxy (row) via iterators") {
        auto m2 = m1;
        auto m3 = m2({0,1,1}, {0,4,1});
        auto vec = std::vector<int> {1,2,3,4};
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m3.begin(), m3.end(), vec.begin()));
    }

    SECTION("Write access to MatrixProxy (column) via iterators") {
        auto m2 = m1;
        auto m3 = m2({0,1,1}, {0,4,1});
        auto vec = std::vector<int> {101,102,103,104};
        std::copy(vec.begin(), vec.end(), m3.begin());
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        auto vec2 = std::vector<int> {101,102,103,104,5,6,7,8,9,10,11,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));
    }

    SECTION("Read individual MatrixProxy (sub-Matrix) items") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        REQUIRE(m3(0,0) == 6);
        REQUIRE(m3(0,1) == 7);
        REQUIRE(m3(1,0) == 10);
        REQUIRE(m3(1,1) == 11);
    }

    SECTION("Write individual MatrixProxy (sub-Matrix) items") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        m3(0,0) = 106;
        REQUIRE(m3(0,0) == 106);

        m3(0,1) = 107;
        REQUIRE(m3(0,1) == 107);

        m3(1,0) = 1010;
        REQUIRE(m3(1,0) == 1010);

        m3(1,1) = 1011;
        REQUIRE(m3(1,1) == 1011);

        auto vec = std::vector<int> {1,2,3,4,5,106,107,8,9,1010,1011,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));
    }

    SECTION("Read access to MatrixProxy (sub-Matrix) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});
        auto vec = std::vector<int> {6,7,10,11};
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m3.begin(), m3.end(), vec.begin()));
    }

    SECTION("Write access to MatrixProxy (sub-Matrix) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});
        auto vec = std::vector<int> {106,107,1010,1011};
        std::copy(vec.begin(), vec.end(), m3.begin());
        REQUIRE(std::equal(m3.begin(), m3.end(), vec.begin()));

        auto vec2 = std::vector<int> {1,2,3,4,5,106,107,8,9,1010,1011,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));
    }

    SECTION("Read individual MatrixProxy (sub-Matrix column) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({0,3,1}, {0,1,1});
        REQUIRE(m4(0,0) == 6);
        REQUIRE(m4(1,0) == 10);
        REQUIRE(m4(2,0) == 14);

        m4 = m3({0,3,1}, {1,1,1});
        REQUIRE(m4(0,0) == 7);
        REQUIRE(m4(1,0) == 11);
        REQUIRE(m4(2,0) == 15);

        m4 = m3({0,3,1}, {2,1,1});
        REQUIRE(m4(0,0) == 8);
        REQUIRE(m4(1,0) == 12);
        REQUIRE(m4(2,0) == 16);
    }

    SECTION("Write individual MatrixProxy (sub-Matrix column) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({0,3,1}, {0,1,1});
        m4(0,0) = 106;
        m4(1,0) = 1010;
        m4(2,0) = 1014;

        REQUIRE(m4(0,0) == 106);
        REQUIRE(m4(1,0) == 1010);
        REQUIRE(m4(2,0) == 1014);
        auto vec = std::vector<int> {1,2,3,4,5,106,7,8,9,1010,11,12,13,1014,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));

        m4 = m3({0,3,1}, {1,1,1});
        m4(0,0) = 107;
        m4(1,0) = 1011;
        m4(2,0) = 1015;

        REQUIRE(m4(0,0) == 107);
        REQUIRE(m4(1,0) == 1011);
        REQUIRE(m4(2,0) == 1015);
        auto vec2 = std::vector<int> {1,2,3,4,5,106,107,8,9,1010,1011,12,13,1014,1015,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));

        m4 = m3({0,3,1}, {2,1,1});
        m4(0,0) = 108;
        m4(1,0) = 1012;
        m4(2,0) = 1016;

        REQUIRE(m4(0,0) == 108);
        REQUIRE(m4(1,0) == 1012);
        REQUIRE(m4(2,0) == 1016);
        auto vec3 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec3.begin()));
    }

    SECTION("Read access to MatrixProxy (sub-Matrix column) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        auto m4 = m3({0,3,1}, {0,1,1});
        auto vec = std::vector<int> {6,10,14};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec.begin()));

        m4 = m3({0,3,1}, {1,1,1});
        auto vec2 = std::vector<int> {7,11,15};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec2.begin()));

        std::reverse(vec2.begin(), vec2.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec2.begin()));

        m4 = m3({0,3,1}, {2,1,1});
        auto vec3 = std::vector<int> {8,12,16};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec3.begin()));

        std::reverse(vec3.begin(), vec3.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec3.begin()));
    }

    SECTION("Write access to MatrixProxy (sub-Matrix column) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        auto m4 = m3({0,3,1}, {0,1,1});
        auto vec = std::vector<int> {106,1010,1014};
        std::copy(vec.begin(), vec.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        auto vec2 = std::vector<int> {1,2,3,4,5,106,7,8,9,1010,11,12,13,1014,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));

        m4 = m3({0,3,1}, {1,1,1});
        auto vec3 = std::vector<int> {107,1011,1015};
        std::copy(vec3.begin(), vec3.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec3.begin()));

        auto vec4 = std::vector<int> {1,2,3,4,5,106,107,8,9,1010,1011,12,13,1014,1015,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec4.begin()));

        m4 = m3({0,3,1}, {2,1,1});
        auto vec5 = std::vector<int> {108,1012,1016};
        std::copy(vec5.begin(), vec5.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec5.begin()));

        auto vec6 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec6.begin()));
    }

    SECTION("Read individual MatrixProxy (sub-Matrix row) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({0,1,1}, {0,3,1});
        REQUIRE(m4(0,0) == 6);
        REQUIRE(m4(0,1) == 7);
        REQUIRE(m4(0,2) == 8);

        m4 = m3({1,1,1}, {0,3,1});
        REQUIRE(m4(0,0) == 10);
        REQUIRE(m4(0,1) == 11);
        REQUIRE(m4(0,2) == 12);

        m4 = m3({2,1,1}, {0,3,1});
        REQUIRE(m4(0,0) == 14);
        REQUIRE(m4(0,1) == 15);
        REQUIRE(m4(0,2) == 16);
    }

    SECTION("Write individual MatrixProxy (sub-Matrix column) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({0,1,1}, {0,3,1});
        m4(0,0) = 106;
        m4(0,1) = 107;
        m4(0,2) = 108;

        REQUIRE(m4(0,0) == 106);
        REQUIRE(m4(0,1) == 107);
        REQUIRE(m4(0,2) == 108);
        auto vec = std::vector<int> {1,2,3,4,5,106,107,108,9,10,11,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));

        m4 = m3({1,1,1}, {0,3,1});
        m4(0,0) = 1010;
        m4(0,1) = 1011;
        m4(0,2) = 1012;

        REQUIRE(m4(0,0) == 1010);
        REQUIRE(m4(0,1) == 1011);
        REQUIRE(m4(0,2) == 1012);
        auto vec2 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));

        m4 = m3({2,1,1}, {0,3,1});
        m4(0,0) = 1014;
        m4(0,1) = 1015;
        m4(0,2) = 1016;

        REQUIRE(m4(0,0) == 1014);
        REQUIRE(m4(0,1) == 1015);
        REQUIRE(m4(0,2) == 1016);
        auto vec3 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec3.begin()));
    }

    SECTION("Read access to MatrixProxy (sub-Matrix row) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        auto m4 = m3({0,1,1}, {0,3,1});
        auto vec = std::vector<int> {6,7,8};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec.begin()));

        m4 = m3({1,1,1}, {0,3,1});
        auto vec2 = std::vector<int> {10,11,12};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec2.begin()));

        std::reverse(vec2.begin(), vec2.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec2.begin()));

        m4 = m3({2,1,1}, {0,3,1});
        auto vec3 = std::vector<int> {14,15,16};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec3.begin()));

        std::reverse(vec3.begin(), vec3.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec3.begin()));
    }

    SECTION("Write access to MatrixProxy (sub-Matrix row) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,2,1}, {1,2,1});

        auto m4 = m3({0,1,1}, {0,3,1});
        auto vec = std::vector<int> {106,107,108};
        std::copy(vec.begin(), vec.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        auto vec2 = std::vector<int> {1,2,3,4,5,106,107,108,9,10,11,12,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));

        m4 = m3({1,1,1}, {0,3,1});
        auto vec3 = std::vector<int> {1010,1011,1012};
        std::copy(vec3.begin(), vec3.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec3.begin()));

        auto vec4 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,14,15,16};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec4.begin()));

        m4 = m3({2,1,1}, {0,3,1});
        auto vec5 = std::vector<int> {1014,1015,1016};
        std::copy(vec5.begin(), vec5.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec5.begin()));

        auto vec6 = std::vector<int> {1,2,3,4,5,106,107,108,9,1010,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec6.begin()));
    }

    SECTION("Read individual MatrixProxy (sub-sub-Matrix) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({1,2,1}, {1,2,1});
        REQUIRE(m4(0,0) == 11);
        REQUIRE(m4(0,1) == 12);
        REQUIRE(m4(1,0) == 15);
        REQUIRE(m4(1,1) == 16);

        m4 = m3({0,2,2}, {0,2,2});
        REQUIRE(m4(0,0) == 6);
        REQUIRE(m4(0,1) == 8);
        REQUIRE(m4(1,0) == 14);
        REQUIRE(m4(1,1) == 16);
    }

    SECTION("Write individual MatrixProxy (sub-sub-Matrix) items") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({1,2,1}, {1,2,1});
        m4(0,0) = 1011;
        m4(0,1) = 1012;
        m4(1,0) = 1015;
        m4(1,1) = 1016;

        REQUIRE(m4(0,0) == 1011);
        REQUIRE(m4(0,1) == 1012);
        REQUIRE(m4(1,0) == 1015);
        REQUIRE(m4(1,1) == 1016);
        auto vec = std::vector<int> {1,2,3,4,5,6,7,8,9,10,1011,1012,13,14,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec.begin()));

        m4 = m3({0,2,2}, {0,2,2});
        m4(0,0) = 106;
        m4(0,1) = 108;
        m4(1,0) = 1014;
        m4(1,1) = 1016;

        REQUIRE(m4(0,0) == 106);
        REQUIRE(m4(0,1) == 108);
        REQUIRE(m4(1,0) == 1014);
        REQUIRE(m4(1,1) == 1016);
        auto vec2 = std::vector<int> {1,2,3,4,5,106,7,108,9,10,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));
    }

    SECTION("Read access to MatrixProxy (sub-sub-Matrix) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({1,2,1}, {1,2,1});
        auto vec = std::vector<int> {11,12,15,16};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        std::reverse(vec.begin(), vec.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec.begin()));

        m4 = m3({0,2,2}, {0,2,2});
        auto vec2 = std::vector<int> {6,8,14,16};
        REQUIRE(std::equal(m4.begin(), m4.end(), vec2.begin()));

        std::reverse(vec2.begin(), vec2.end());
        REQUIRE_FALSE(std::equal(m4.begin(), m4.end(), vec2.begin()));
    }

    SECTION("Write access to MatrixProxy (sub-sub-Matrix) via iterators") {
        auto m2 = m1;
        auto m3 = m2({1,3,1}, {1,3,1});

        auto m4 = m3({1,2,1}, {1,2,1});
        auto vec = std::vector<int> {1011,1012,1015,1016};
        std::copy(vec.begin(), vec.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec.begin()));

        auto vec2 = std::vector<int> {1,2,3,4,5,6,7,8,9,10,1011,1012,13,14,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec2.begin()));

        m4 = m3({0,2,2}, {0,2,2});
        auto vec3 = std::vector<int> {106,108,1014, 1016};
        std::copy(vec3.begin(), vec3.end(), m4.begin());
        REQUIRE(std::equal(m4.begin(), m4.end(), vec3.begin()));

        auto vec4 = std::vector<int> {1,2,3,4,5,106,7,108,9,10,1011,1012,13,1014,1015,1016};
        REQUIRE(std::equal(m2.begin(), m2.end(), vec4.begin()));
    }

}
