set(ASL_SOURCES
        ${ASL_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/gjk.c
        ${CMAKE_CURRENT_SOURCE_DIR}/bvh.c
        ${CMAKE_CURRENT_SOURCE_DIR}/bvh_f16c.c
        ${CMAKE_CURRENT_SOURCE_DIR}/moeller.c
        PARENT_SCOPE
)