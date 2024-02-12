//
// Created by kenne on 11/02/2024.
//

#pragma once

namespace nxx::roots {

  template<std::integral ITER_T, typename... ARGS>
  using IterData = std::tuple<ITER_T, ARGS...>;

}
