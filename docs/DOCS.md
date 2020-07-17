<h1 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/logo.svg?sanitize=true"  alt="Non-Repetitive Parts Calculator" width="418" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Background">Background</a> •
  <a href="#Finder-Mode">Finder Mode</a> •
  <a href="#Maker-Mode">Maker Mode</a> •
  <a href="#NRP-Calculator-in-Action">NRP Calculator in Action</a> •
  <a href="../README.md">README</a>
</p>

### `ShareDB` API Documentation
---
**\_\_init__(self, path, reset=False, serial='msgpack', compress=False, readers=100, buffer_size=10\*\*5, map_size=10\*\*9)**

`ShareDB` **constructor**.

| argument | type | description | default |
|--|--|--|--|
| `path` | `string` | a/path/to/a/directory/to/persist/the/data |  -- |
| `reset` | `boolean` | if `True` - delete and recreate path following subsequent parameters | `False` |
| `serial` | `string` | must be either `'msgpack'` or `'pickle'` | `'msgpack'` |
| `compress` | `string` | if `True` - will compress the values using `zlib` | `False` |
| `readers` | `integer` | max no. of processes that may read data in parallel | `100` |
| `buffer_size` | `integer` | max no. of commits after which a sync is triggered | `100,000` |
| `map_size` | `integer` | max amount of bytes to allocate for storage | `10**9` (1 GB) |

**_Returns_**: `self` to `ShareDB` object.

```python
>>> from ShareDB import ShareDB
>>> myDB = ShareDB(
    path='./test.ShareDB',
    reset=True,
    readers=10,
    buffer_size=100,
    map_size=10**5)
>>> myDB.ALIVE
True
```
---