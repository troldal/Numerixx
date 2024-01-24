//
// Created by kenne on 08/01/2024.
//

#pragma once

#include "_external.hpp"

namespace nxx::traits
{
    /**
     * @brief A trait class to deduce the value type of a given container.
     *
     * @details This structure template is used to determine the value type of a container.
     *          It defaults to using the `value_type` member of the container, if available.
     *          Specializations can be provided for containers that do not follow this standard.
     *
     * @tparam Container The container type for which the value type is being determined.
     * @tparam <unnamed> An optional, unnamed template parameter for SFINAE.
     */
    template<typename Container, typename = std::void_t<>>
    struct ContainerValueType
    {
        using type = typename Container::value_type;
    };

    /**
     * @brief Specialization of ContainerValueType for blaze::DynamicVector.
     *
     * @details This specialization handles the value type deduction for blaze::DynamicVector.
     *          It directly uses the template parameter `T` as the value type.
     *
     * @tparam T The type of elements in the blaze::DynamicVector.
     * @tparam SO The storage order of the blaze::DynamicVector.
     */
    template<typename T, bool SO>
    struct ContainerValueType<blaze::DynamicVector<T, SO>, std::void_t<>> {
        using type = T;
    };

    /**
     * @brief Alias template for easier access to the value type of a container.
     *
     * @details Provides a convenient way to retrieve the value type of a container
     *          using ContainerValueType, without having to access its `type` member directly.
     *
     * @tparam Container The container type for which the value type is deduced.
     */
    template<typename Container>
    using ContainerValueType_t = typename ContainerValueType<Container>::type;
}

