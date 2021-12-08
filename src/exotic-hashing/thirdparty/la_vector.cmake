include(FetchContent)

set(LA_VECTOR_LIBRARY la_vector)
FetchContent_Declare(
  ${LA_VECTOR_LIBRARY}
  GIT_REPOSITORY https://github.com/DominikHorn/la_vector.git
  GIT_TAG cd0409e
  )

FetchContent_MakeAvailable(${LA_VECTOR_LIBRARY})
