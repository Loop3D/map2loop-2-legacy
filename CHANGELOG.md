# Changelog

## [1.3.3](https://github.com/Loop3D/map2loop-2/compare/1.3.2...1.3.3) (2022-12-04)


### Bug Fixes

* version bump ([d9bdd34](https://github.com/Loop3D/map2loop-2/commit/d9bdd343461f1bf278da4e7665aca18ebf5a8422))

## [1.3.2](https://github.com/Loop3D/map2loop-2/compare/1.3.1...1.3.2) (2022-11-30)


### Bug Fixes

* Add clipping of local file dtm ([3359927](https://github.com/Loop3D/map2loop-2/commit/3359927b8ef7f9bc78527cb0e6f70982deca581f))
* Fix strike conversion logic, Replace drift_prefix with ignore_codes option, remove opentopography dtm ([97d26ed](https://github.com/Loop3D/map2loop-2/commit/97d26ed0f6f6d5545a4c3c3ad2c47dee3bfe7c9a))
* Move sjoin out of loop for efficiency ([fe59514](https://github.com/Loop3D/map2loop-2/commit/fe59514b2737aaed751b334bd2a62bd47b8d7f5b))
* version bump ([e9200d1](https://github.com/Loop3D/map2loop-2/commit/e9200d18337835546aea202f0ac62ec705253ad2))

## [1.3.1](https://github.com/Loop3D/map2loop-2/compare/1.3.0...1.3.1) (2022-11-08)


### Bug Fixes

* version bump for mambaforge build ([31eabb3](https://github.com/Loop3D/map2loop-2/commit/31eabb37eab016b30d918d8793f159769fa4313c))

## [1.3.0](https://github.com/Loop3D/map2loop-2/compare/v1.2.10...1.3.0) (2022-11-01)


### Features

* Merged classes into master, with version bump ([4cc796e](https://github.com/Loop3D/map2loop-2/commit/4cc796e38abcc4d7be4a641b64da54b17e637eed))


### Bug Fixes

* add crs to temporary geodataframe to solve crs mismatch ([4d7023a](https://github.com/Loop3D/map2loop-2/commit/4d7023aa45265b5b729d5cea573fe955103ba6e3))
* Add supergroup to stratigraphic column, fix thickness values in projectfile output ([8387655](https://github.com/Loop3D/map2loop-2/commit/83876550293d0143d6dac0914d67807e93343bad))
* Move crs assignment to after geology has been specified ([e74e9e9](https://github.com/Loop3D/map2loop-2/commit/e74e9e91918571b227ed368e8ef21b21f3906fe7))
* rework of topology file structure and graph usage ([a5419e8](https://github.com/Loop3D/map2loop-2/commit/a5419e8d8f70a5dbd41db4fd842a0127b52b89aa))
* Without explicit string type replacement silently fails ([451a876](https://github.com/Loop3D/map2loop-2/commit/451a876b9c3549e49d3f3f9691d87c2ee595ddd0))

## [1.2.10](https://github.com/Loop3D/map2loop-2/compare/v1.2.9...v1.2.10) (2022-08-20)


### Bug Fixes

* Remove pypi dependencies while gdal and rasterio can't install easily ([17c9076](https://github.com/Loop3D/map2loop-2/commit/17c9076d11e104e24ab79aa777c69e951f0c058c))

### [1.2.9](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.8...v1.2.9) (2022-08-18)


### Bug Fixes

* Reset tags to not include v and attempt to get pip install working ([2d0924b](https://www.github.com/Loop3D/map2loop-2/commit/2d0924b4127e7e659e99a2c3258b5cb983e28bd0))

### [1.2.8](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.7...v1.2.8) (2022-08-03)


### Bug Fixes

* change rasterio dependency to enable python 3.7 to work ([2fb5b80](https://www.github.com/Loop3D/map2loop-2/commit/2fb5b80ad9c86afae9b3d79fef2acfceea4e0d43))

### [1.2.7](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.6...v1.2.7) (2022-08-02)


### Bug Fixes

* ordering ([8901d03](https://www.github.com/Loop3D/map2loop-2/commit/8901d0306ba9c1659869662449af76ae521aa324))

### [1.2.6](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.5...v1.2.6) (2022-08-02)


### Bug Fixes

* remove osx conda build ([7209704](https://www.github.com/Loop3D/map2loop-2/commit/72097049f479714c9c3e73941bb978e4523977b4))

### [1.2.5](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.4...v1.2.5) (2022-08-02)


### Bug Fixes

* Add python version to continuous integration ([6835882](https://www.github.com/Loop3D/map2loop-2/commit/6835882ab1f3723e016ed7fd4bee71419be996fc))

### [1.2.4](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.3...v1.2.4) (2022-08-02)


### Bug Fixes

* Temporary fix to avoid rasterio/gdal conflicts ([23e5f6a](https://www.github.com/Loop3D/map2loop-2/commit/23e5f6a9fce11281add61161fb0d6637b37d1da9))

### [1.2.3](https://www.github.com/Loop3D/map2loop-2/compare/v1.2.2...v1.2.3) (2022-07-27)


### Bug Fixes

* added function to check if the filenames exist ([f7b6475](https://www.github.com/Loop3D/map2loop-2/commit/f7b6475333a74895644a3f6b19dcb3bd750f5fc0))
* black format ([314ee87](https://www.github.com/Loop3D/map2loop-2/commit/314ee87f5f8a266ac2ab724f0b3faff2bb89e20d))
* casting numeric columns as float ([2755546](https://www.github.com/Loop3D/map2loop-2/commit/275554669d70c5319d81fd58c1afd2497d02685c))
* normalising fault trace vector before comparing to dip direction ([bda0ddd](https://www.github.com/Loop3D/map2loop-2/commit/bda0dddfd6d4dee0cc6b347d62ac08da540671db))
* removing type hint for python 3.8 ([587441f](https://www.github.com/Loop3D/map2loop-2/commit/587441fb87b7fdf34ef6abe66cc1b771cf0d94c8))


### Documentation

* adding templates for issues ([2874b97](https://www.github.com/Loop3D/map2loop-2/commit/2874b97feb883078a5e7af7a8a86fd78ec76bada))
