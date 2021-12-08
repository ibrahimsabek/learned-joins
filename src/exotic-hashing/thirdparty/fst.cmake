include(FetchContent)

set(FST_LIBRARY fst)
FetchContent_Declare(
  ${FST_LIBRARY}
  GIT_REPOSITORY https://github.com/DominikHorn/FST
  GIT_TAG c21fbaa
)

FetchContent_MakeAvailable(${FST_LIBRARY})
